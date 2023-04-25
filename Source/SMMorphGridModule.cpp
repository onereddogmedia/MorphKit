// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMMorphGridModule.h"
#include "SMLeakDebugger.h"
#include "SMLiveDecoder.h"
#include "SMMath.h"
#include "SMMorphGrid.h"
#include "SMMorphPlanVoice.h"
#include "SMMorphUtils.h"

#include <assert.h>

using namespace SpectMorph;

using std::max;
using std::min;
using std::sort;
using std::string;
using std::vector;

static LeakDebugger leak_debugger("SpectMorph::MorphGridModule");

MorphGridModule::MorphGridModule(MorphPlanVoice* voice) : MorphOperatorModule(voice) {
    leak_debugger.add(this);

    my_source.module = this;

    audio.fundamental_freq = 440;
    audio.mix_freq = 48000;
    audio.frame_size_ms = 1;
    audio.frame_step_ms = 1;
    audio.zeropad = 4;
    audio.loop_type = Audio::LOOP_NONE;
}

MorphGridModule::~MorphGridModule() {
    leak_debugger.del(this);
}

void MorphGridModule::set_config(const MorphOperatorConfig* op_cfg) {
    cfg = dynamic_cast<const MorphGrid::Config*>(op_cfg);
    g_return_if_fail(cfg != nullptr);

    input_node.resize(cfg->width);

    for (unsigned int x = 0; x < cfg->width; x++) {
        input_node[(size_t)x].resize(cfg->height);
        for (unsigned int y = 0; y < cfg->height; y++) {
            const MorphGridNode& node = cfg->input_node[x][y];

            if (node.path != "") {
                input_node[x][y].source.set_wav_set(cfg->wav_set_repo, node.path);
                input_node[x][y].has_source = true;
            } else {
                input_node[x][y].has_source = false;
            }
        }
    }

    clear_dependencies();
}

void MorphGridModule::MySource::prepareToPlay(float mix_freq) {
    for (unsigned int x = 0; x < module->cfg->width; x++) {
        for (unsigned int y = 0; y < module->cfg->height; y++) {
            InputNode& node = module->input_node[x][y];
            if (node.has_source) {
                node.source.prepareToPlay(mix_freq);
            }
        }
    }
}

void MorphGridModule::MySource::retrigger(int channel, float freq, int midi_velocity) {
    for (unsigned int x = 0; x < module->cfg->width; x++) {
        for (unsigned int y = 0; y < module->cfg->height; y++) {
            InputNode& node = module->input_node[x][y];
            if (node.has_source) {
                node.source.retrigger(channel, freq, midi_velocity);
            }
        }
    }
}

Audio* MorphGridModule::MySource::audio() {
    return &module->audio;
}

static bool get_normalized_block(MorphGridModule::InputNode& input_node, size_t index, AudioBlock& out_audio_block) {
    LiveDecoderSource* source = nullptr;

    if (input_node.has_source) {
        source = &input_node.source;
    }
    const double time_ms = index; // 1ms frame step

    return MorphUtils::get_normalized_block(source, time_ms, out_audio_block);
}

namespace {
struct MagData {
    enum { BLOCK_LEFT = 0, BLOCK_RIGHT = 1 } block;
    size_t index;
    uint16_t mag;
};

static bool md_cmp(const MagData& m1, const MagData& m2) {
    return m1.mag > m2.mag; // sort with biggest magnitude first
}

static void interp_mag_one(double interp, uint16_t* left, uint16_t* right) {
    const uint16_t lmag_idb = max<uint16_t>(left ? *left : 0, SM_IDB_CONST_M96);
    const uint16_t rmag_idb = max<uint16_t>(right ? *right : 0, SM_IDB_CONST_M96);

    const uint16_t mag_idb = (uint16_t)sm_round_positive((1 - interp) * lmag_idb + interp * rmag_idb);

    if (left)
        *left = mag_idb;
    if (right)
        *right = mag_idb;
}

static void morph_scale(AudioBlock& out_block, const AudioBlock& in_block, double factor) {
    const int ddb = sm_factor2delta_idb(factor);

    out_block = in_block;
    for (size_t i = 0; i < out_block.noise.size(); i++)
        out_block.noise[i] = (uint16_t)sm_bound<int>(0, out_block.noise[i] + ddb, 65535);

    for (size_t i = 0; i < out_block.freqs.size(); i++)
        interp_mag_one(factor, nullptr, &out_block.mags[i]);
}

} // namespace

static bool morph(AudioBlock& out_block, bool have_left, const AudioBlock& left_block, bool have_right,
                  const AudioBlock& right_block, double morphing) {
    const double interp = (morphing + 1) / 2; /* examples => 0: only left; 0.5 both equally; 1: only right */

    if (!have_left && !have_right) // nothing + nothing = nothing
        return false;

    if (!have_left) // nothing + interp * right = interp * right
    {
        morph_scale(out_block, right_block, interp);
        return true;
    }
    if (!have_right) // (1 - interp) * left + nothing = (1 - interp) * left
    {
        morph_scale(out_block, left_block, 1 - interp);
        return true;
    }

    // clear result block
    out_block.freqs.clear();
    out_block.mags.clear();

    // FIXME: lpc stuff
    MagData mds[AudioBlock::block_size * 2];
    size_t mds_size = 0;
    for (size_t i = 0; i < left_block.freqs.size(); i++) {
        MagData& md = mds[mds_size];

        md.block = MagData::BLOCK_LEFT;
        md.index = i;
        md.mag = left_block.mags[i];
        mds_size++;
    }
    for (size_t i = 0; i < right_block.freqs.size(); i++) {
        MagData& md = mds[mds_size];

        md.block = MagData::BLOCK_RIGHT;
        md.index = i;
        md.mag = right_block.mags[i];
        mds_size++;
    }
    sort(mds, mds + mds_size, md_cmp);

    size_t left_freqs_size = left_block.freqs.size();
    size_t right_freqs_size = right_block.freqs.size();

    MorphUtils::FreqState left_freqs[AudioBlock::block_size];
    MorphUtils::FreqState right_freqs[AudioBlock::block_size];

    init_freq_state(left_block.freqs, left_freqs);
    init_freq_state(right_block.freqs, right_freqs);

    for (size_t m = 0; m < mds_size; m++) {
        size_t i = 0, j = 0;
        bool match = false;
        if (mds[m].block == MagData::BLOCK_LEFT) {
            i = mds[m].index;

            if (!left_freqs[i].used)
                match = MorphUtils::find_match(left_freqs[i].freq_f, right_freqs, right_freqs_size, &j);
        } else // (mds[m].block == MagData::BLOCK_RIGHT)
        {
            j = mds[m].index;
            if (!right_freqs[j].used)
                match = MorphUtils::find_match(right_freqs[j].freq_f, left_freqs, left_freqs_size, &i);
        }
        if (match) {
            /* prefer frequency of louder partial:
             *
             * if the magnitudes are similar, mfact will be close to 1, and freq will become approx.
             *
             *   freq = (1 - interp) * lfreq + interp * rfreq
             *
             * if the magnitudes are very different, mfact will be close to 0, and freq will become
             *
             *   freq ~= lfreq         // if left partial is louder
             *   freq ~= rfreq         // if right partial is louder
             */
            const double lfreq = left_block.freqs[i];
            const double rfreq = right_block.freqs[j];
            double freq;

            if (left_block.mags[i] > right_block.mags[j]) {
                const double mfact = right_block.mags_f(j) / left_block.mags_f(i);

                freq = lfreq + mfact * interp * (rfreq - lfreq);
            } else {
                const double mfact = left_block.mags_f(i) / right_block.mags_f(j);

                freq = rfreq + mfact * (1 - interp) * (lfreq - rfreq);
            }
            // FIXME: lpc
            // FIXME: non-db

            const uint16_t lmag_idb = max(left_block.mags[i], SM_IDB_CONST_M96);
            const uint16_t rmag_idb = max(right_block.mags[j], SM_IDB_CONST_M96);
            const uint16_t mag_idb = (uint16_t)sm_round_positive((1 - interp) * lmag_idb + interp * rmag_idb);

            out_block.freqs.push_back((unsigned short)freq);
            out_block.mags.push_back(mag_idb);

            left_freqs[i].used = 1;
            right_freqs[j].used = 1;
        }
    }
    for (size_t i = 0; i < left_freqs_size; i++) {
        if (!left_freqs[i].used) {
            out_block.freqs.push_back(left_block.freqs[i]);
            out_block.mags.push_back(left_block.mags[i]);

            interp_mag_one(interp, &out_block.mags.back(), nullptr);
        }
    }
    for (size_t i = 0; i < right_freqs_size; i++) {
        if (!right_freqs[i].used) {
            out_block.freqs.push_back(right_block.freqs[i]);
            out_block.mags.push_back(right_block.mags[i]);

            interp_mag_one(interp, nullptr, &out_block.mags.back());
        }
    }
    out_block.noise.clear();
    for (size_t i = 0; i < left_block.noise.size(); i++)
        out_block.noise.push_back(
            sm_factor2idb((1 - interp) * left_block.noise_f(i) + interp * right_block.noise_f(i)));

    out_block.sort_freqs();
    return true;
}

namespace {

struct LocalMorphParams {
    unsigned int start;
    unsigned int end;
    double morphing;
};

static LocalMorphParams global_to_local_params(double global_morphing, unsigned int node_count) {
    LocalMorphParams result;

    /* interp: range for node_count=3: 0 ... 2.0 */
    const double interp = (global_morphing + 1) / 2 * (node_count - 1);

    // find the two adjecant nodes (double -> integer position)
    result.start = (uint16_t)sm_bound<int>(0, (int)interp, (int)(node_count - 1));
    result.end = (uint16_t)sm_bound<int>(0, (int)result.start + 1, (int)(node_count - 1));

    if (approximatelyEqual(interp, (double)result.end)) {
        result.morphing = result.end;
        return result;
    }
    const double interp_frac = sm_bound(0.0, interp - result.start, 1.0); /* position between adjecant nodes */
    result.morphing = interp_frac * 2 - 1; /* normalize fractional part to range -1.0 ... 1.0 */
    return result;
}

static void apply_delta_db(AudioBlock& block, double delta_db) {
    double factor;
    if (delta_db <= -12) {
        factor = -96;
    } else {
        factor = db_to_factor(delta_db);
    }
    const int ddb = sm_factor2delta_idb(factor);

    // apply delta db volume to partials & noise
    for (size_t i = 0; i < block.mags.size(); i++)
        block.mags[i] = (uint16_t)sm_bound<int>(0, block.mags[i] + ddb, 65535);

    for (size_t i = 0; i < block.noise.size(); i++)
        block.noise[i] = (uint16_t)sm_bound<int>(0, block.noise[i] + ddb, 65535);
}

} // namespace

static double get_morphing(MorphGrid::ControlType type, double gui_value, MorphPlanVoice* voice) {
    return voice->control_input(gui_value, type);
}

AudioBlock* MorphGridModule::MySource::audio_block(size_t index) {
    const double x_morphing =
        get_morphing(module->cfg->x_control_type, module->cfg->x_morphing, module->morph_plan_voice);
    const double y_morphing =
        get_morphing(module->cfg->y_control_type, module->cfg->y_morphing, module->morph_plan_voice);

    double a_db = get_morphing(module->cfg->node_a_db_control_type, 0, module->morph_plan_voice);
    double b_db = get_morphing(module->cfg->node_b_db_control_type, 0, module->morph_plan_voice);
    double c_db = get_morphing(module->cfg->node_c_db_control_type, 0, module->morph_plan_voice);
    double d_db = get_morphing(module->cfg->node_d_db_control_type, 0, module->morph_plan_voice);

    const LocalMorphParams x_morph_params = global_to_local_params(x_morphing, module->cfg->width);
    const LocalMorphParams y_morph_params = global_to_local_params(y_morphing, module->cfg->height);

    double sum_a = 0;
    double sum_b = 0;
    double sum_c = 0;
    double sum_d = 0;

    if (module->cfg->width == 1 && module->cfg->height == 1) {
        /*
         *  A (UI D)
         */
        InputNode& node_a = module->input_node[0][y_morph_params.start];
        get_normalized_block(node_a, index, audio_block_a);
        apply_delta_db(audio_block_a, d_db);
        module->audio_block = audio_block_a;

        for (size_t i = 0; i < audio_block_a.mags.size(); i++) {
            double s = sm_idb2factor_slow(audio_block_a.mags[i]);
            sum_d += s * s;
        }
    } else if (module->cfg->width == 1) {
        /*
         *  A (UI A)
         *  |
         *  |
         *  B (UI D)
         */
        InputNode& node_a = module->input_node[0][y_morph_params.start];
        InputNode& node_b = module->input_node[0][y_morph_params.end];

        bool have_a = get_normalized_block(node_a, index, audio_block_a);
        bool have_b = get_normalized_block(node_b, index, audio_block_b);

        apply_delta_db(audio_block_a, b_db); // D
        apply_delta_db(audio_block_b, a_db); // A

        morph(module->audio_block, have_a, audio_block_a, have_b, audio_block_b, y_morph_params.morphing);

        for (size_t i = 0; i < audio_block_a.mags.size(); i++) {
            double s = sm_idb2factor_slow(audio_block_a.mags[i]);
            sum_a += s * s;
        }
        for (size_t i = 0; i < audio_block_b.mags.size(); i++) {
            double s = sm_idb2factor_slow(audio_block_b.mags[i]);
            sum_d += s * s;
        }
    } else {
        /*                UI
         *  A ---- B   A ---- B
         *  |      |   |      |
         *  |      |   |      |
         *  C ---- D   D ---- C
         */
        InputNode& node_a = module->input_node[x_morph_params.start][y_morph_params.start];
        InputNode& node_b = module->input_node[x_morph_params.end][y_morph_params.start];
        InputNode& node_c = module->input_node[x_morph_params.start][y_morph_params.end];
        InputNode& node_d = module->input_node[x_morph_params.end][y_morph_params.end];

        bool have_a = get_normalized_block(node_a, index, audio_block_a);
        bool have_b = get_normalized_block(node_b, index, audio_block_b);
        bool have_c = get_normalized_block(node_c, index, audio_block_c);
        bool have_d = get_normalized_block(node_d, index, audio_block_d);

        apply_delta_db(audio_block_a, b_db); // D    (UI)
        apply_delta_db(audio_block_b, d_db); // C
        apply_delta_db(audio_block_c, a_db); // A
        apply_delta_db(audio_block_d, c_db); // B

        bool have_ab = morph(audio_block_ab, have_a, audio_block_a, have_b, audio_block_b, x_morph_params.morphing);
        bool have_cd = morph(audio_block_cd, have_c, audio_block_c, have_d, audio_block_d, x_morph_params.morphing);
        morph(module->audio_block, have_ab, audio_block_ab, have_cd, audio_block_cd, y_morph_params.morphing);

        for (size_t i = 0; i < audio_block_a.mags.size(); i++) {
            double s = sm_idb2factor_slow(audio_block_a.mags[i]);
            sum_a += s * s;
        }
        for (size_t i = 0; i < audio_block_b.mags.size(); i++) {
            double s = sm_idb2factor_slow(audio_block_b.mags[i]);
            sum_b += s * s;
        }
        for (size_t i = 0; i < audio_block_c.mags.size(); i++) {
            double s = sm_idb2factor_slow(audio_block_c.mags[i]);
            sum_c += s * s;
        }
        for (size_t i = 0; i < audio_block_d.mags.size(); i++) {
            double s = sm_idb2factor_slow(audio_block_d.mags[i]);
            sum_d += s * s;
        }
    }

    double rms_a = sqrt(sum_a / audio_block_a.mags.size());
    double rms_b = sqrt(sum_b / audio_block_b.mags.size());
    double rms_c = sqrt(sum_c / audio_block_c.mags.size());
    double rms_d = sqrt(sum_d / audio_block_d.mags.size());

    double mag_a, mag_b, mag_c, mag_d;
    double acmix, bdmix;
    calculateVectorMixValues<double>(0, 0, x_morphing, y_morphing, mag_a, mag_b, mag_c, mag_d, acmix, bdmix, 2, false);

    // round off very small values
    if (mag_a < 0.01)
        mag_a = 0;
    if (mag_b < 0.01)
        mag_b = 0;
    if (mag_c < 0.01)
        mag_c = 0;
    if (mag_d < 0.01)
        mag_d = 0;

    module->morph_plan_voice->set_control_output(0, rms_a * mag_a);
    module->morph_plan_voice->set_control_output(1, rms_b * mag_b);
    module->morph_plan_voice->set_control_output(2, rms_c * mag_c);
    module->morph_plan_voice->set_control_output(3, rms_d * mag_d);

    return &module->audio_block;
}

LiveDecoderSource* MorphGridModule::source() {
    return &my_source;
}
