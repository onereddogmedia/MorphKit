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

    for (int x = 0; x < cfg->width; x++) {
        for (int y = 0; y < cfg->height; y++) {
            const MorphGridNode& node = cfg->input_node[(size_t)x][(size_t)y];

            input_nodes(x, y).mod = morph_plan_voice->module(node.op);

            if (node.wav_set) {
                input_nodes(x, y).source.set_wav_set(node.wav_set);
                input_nodes(x, y).has_source = true;
            } else {
                input_nodes(x, y).has_source = false;
            }
        }
    }
}

void MorphGridModule::MySource::prepareToPlay(float) {
    // pj probably not needed now
}

void MorphGridModule::MySource::retrigger(int channel, float freq, int midi_velocity, bool onset) {
    for (int x = 0; x < module->cfg->width; x++) {
        for (int y = 0; y < module->cfg->height; y++) {
            InputNode& node = module->input_nodes(x, y);

            if (node.mod && node.mod->source()) {
                node.mod->source()->retrigger(channel, freq, midi_velocity, onset);
            }
            if (node.has_source) {
                node.source.retrigger(channel, freq, midi_velocity, onset);
            }
        }
    }
}

Audio* MorphGridModule::MySource::audio() {
    return &module->audio;
}

static bool get_normalized_block(MorphGridModule::InputNode& input_node, size_t index, RTAudioBlock& out_audio_block) {
    LiveDecoderSource* source = nullptr;

    if (input_node.mod) {
        source = input_node.mod->source();
    } else if (input_node.has_source) {
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

static void morph_scale(RTAudioBlock& out_block, const RTAudioBlock& in_block, double factor) {
    const int ddb = sm_factor2delta_idb(factor);

    out_block.assign(in_block);
    for (size_t i = 0; i < out_block.noise.size(); i++)
        out_block.noise[i] = (uint16_t)sm_bound<int>(0, out_block.noise[i] + ddb, 65535);

    for (size_t i = 0; i < out_block.freqs.size(); i++)
        interp_mag_one(factor, nullptr, &out_block.mags[i]);
}

} // namespace

static bool morph(RTAudioBlock& out_block, bool have_left, const RTAudioBlock& left_block, bool have_right,
                  const RTAudioBlock& right_block, double morphing) {
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

    // set out_block capacity
    const size_t max_partials = left_block.freqs.size() + right_block.freqs.size();
    out_block.freqs.set_capacity(max_partials);
    out_block.mags.set_capacity(max_partials);

    // FIXME: lpc stuff
    MagData mds[max_partials + AVOID_ARRAY_UB];
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

    MorphUtils::FreqState left_freqs[left_freqs_size + AVOID_ARRAY_UB];
    MorphUtils::FreqState right_freqs[right_freqs_size + AVOID_ARRAY_UB];

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
    out_block.noise.set_capacity(left_block.noise.size());
    for (size_t i = 0; i < left_block.noise.size(); i++)
        out_block.noise.push_back(
            sm_factor2idb((1 - interp) * left_block.noise_f(i) + interp * right_block.noise_f(i)));

    out_block.sort_freqs();
    return true;
}

namespace {

struct LocalMorphParams {
    int start;
    int end;
    double morphing;
};

static LocalMorphParams global_to_local_params(double global_morphing, int node_count) {
    LocalMorphParams result;

    /* interp: range for node_count=3: 0 ... 2.0 */
    const double interp = (global_morphing + 1) / 2 * (node_count - 1);

    // find the two adjecant nodes (double -> integer position)
    result.start = (uint16_t)sm_bound<int>(0, (int)interp, (node_count - 1));
    result.end = (uint16_t)sm_bound<int>(0, (int)result.start + 1, (node_count - 1));

    if (approximatelyEqual(interp, (double)result.end)) {
        result.morphing = result.end;
        return result;
    }
    const double interp_frac = sm_bound(0.0, interp - result.start, 1.0); /* position between adjecant nodes */
    result.morphing = interp_frac * 2 - 1; /* normalize fractional part to range -1.0 ... 1.0 */
    return result;
}

static void apply_delta_db(RTAudioBlock& block, double delta_db) {
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

bool MorphGridModule::MySource::rt_audio_block(size_t index, RTAudioBlock& out_block) {
    bool status = false;
    const double x_morphing = module->morph_plan_voice->control_input(0, MorphOperator::CONTROL_SIGNAL_1, module);
    const double y_morphing = module->morph_plan_voice->control_input(0, MorphOperator::CONTROL_SIGNAL_2, module);

    double a_db = module->morph_plan_voice->control_input(0, MorphOperator::CONTROL_SIGNAL_3, module);
    double b_db = module->morph_plan_voice->control_input(0, MorphOperator::CONTROL_SIGNAL_4, module);
    double c_db = module->morph_plan_voice->control_input(0, MorphOperator::CONTROL_SIGNAL_5, module);
    double d_db = module->morph_plan_voice->control_input(0, MorphOperator::CONTROL_SIGNAL_6, module);

    const LocalMorphParams x_morph_params = global_to_local_params(x_morphing, module->cfg->width);
    const LocalMorphParams y_morph_params = global_to_local_params(y_morphing, module->cfg->height);

    double sum_a = 0;
    double sum_b = 0;
    double sum_c = 0;
    double sum_d = 0;
    double rms_a = 0;
    double rms_b = 0;
    double rms_c = 0;
    double rms_d = 0;

    RTAudioBlock audio_block_a(module->rt_memory_area());
    RTAudioBlock audio_block_b(module->rt_memory_area());
    RTAudioBlock audio_block_c(module->rt_memory_area());
    RTAudioBlock audio_block_d(module->rt_memory_area());

    if (module->cfg->width == 1 && module->cfg->height == 1) {
        /*
         *  A (UI D)
         */
        InputNode& node_a = module->input_nodes(x_morph_params.start, 0);
        get_normalized_block(node_a, index, audio_block_a);
        apply_delta_db(audio_block_a, b_db);
        morph(out_block, false, audio_block_a, true, audio_block_a, 1.0f);

        for (size_t i = 0; i < audio_block_a.mags.size(); i++) {
            double s = sm_idb2factor_slow(audio_block_a.mags[i]);
            sum_d += s * s;
        }

        rms_d = sqrt(sum_d / audio_block_a.mags.size());

        status = true;
    } else if (module->cfg->width == 1) {
        /*
         *  A (UI A)
         *  |
         *  |
         *  B (UI D)
         */
        InputNode& node_a = module->input_nodes(0, y_morph_params.start);
        InputNode& node_b = module->input_nodes(0, y_morph_params.end);

        bool have_a = get_normalized_block(node_a, index, audio_block_a);
        bool have_b = get_normalized_block(node_b, index, audio_block_b);

        apply_delta_db(audio_block_a, b_db); // D
        apply_delta_db(audio_block_b, a_db); // A

        morph(out_block, have_a, audio_block_a, have_b, audio_block_b, y_morph_params.morphing);
        for (size_t i = 0; i < audio_block_b.mags.size(); i++) {
            double s = sm_idb2factor_slow(audio_block_b.mags[i]);
            sum_a += s * s;
        }
        for (size_t i = 0; i < audio_block_a.mags.size(); i++) {
            double s = sm_idb2factor_slow(audio_block_a.mags[i]);
            sum_d += s * s;
        }

        rms_a = sqrt(sum_a / audio_block_b.mags.size());
        rms_d = sqrt(sum_d / audio_block_a.mags.size());

        status = true;
    } else {
        RTAudioBlock audio_block_ab(module->rt_memory_area());
        RTAudioBlock audio_block_cd(module->rt_memory_area());

        /*                UI
         *  A ---- B   A ---- B
         *  |      |   |      |
         *  |      |   |      |
         *  C ---- D   D ---- C
         */
        InputNode& node_a = module->input_nodes(x_morph_params.start, y_morph_params.start);
        InputNode& node_b = module->input_nodes(x_morph_params.end, y_morph_params.start);
        InputNode& node_c = module->input_nodes(x_morph_params.start, y_morph_params.end);
        InputNode& node_d = module->input_nodes(x_morph_params.end, y_morph_params.end);

        bool have_a = get_normalized_block(node_a, index, audio_block_a);
        bool have_b = get_normalized_block(node_b, index, audio_block_b);
        bool have_c = get_normalized_block(node_c, index, audio_block_c);
        bool have_d = get_normalized_block(node_d, index, audio_block_d);

        apply_delta_db(audio_block_a, b_db); // D    (UI)
        apply_delta_db(audio_block_b, d_db); // A
        apply_delta_db(audio_block_c, a_db); // C
        apply_delta_db(audio_block_d, c_db); // B

        bool have_ab = morph(audio_block_ab, have_a, audio_block_a, have_b, audio_block_b, x_morph_params.morphing);
        bool have_cd = morph(audio_block_cd, have_c, audio_block_c, have_d, audio_block_d, x_morph_params.morphing);
        bool have_abcd = morph(out_block, have_ab, audio_block_ab, have_cd, audio_block_cd, y_morph_params.morphing);

        if (have_abcd) {
            for (size_t i = 0; i < audio_block_d.mags.size(); i++) {
                double s = sm_idb2factor_slow(audio_block_d.mags[i]);
                sum_a += s * s;
            }
            for (size_t i = 0; i < audio_block_c.mags.size(); i++) {
                double s = sm_idb2factor_slow(audio_block_c.mags[i]);
                sum_b += s * s;
            }
            for (size_t i = 0; i < audio_block_a.mags.size(); i++) {
                double s = sm_idb2factor_slow(audio_block_a.mags[i]);
                sum_c += s * s;
            }
            for (size_t i = 0; i < audio_block_b.mags.size(); i++) {
                double s = sm_idb2factor_slow(audio_block_b.mags[i]);
                sum_d += s * s;
            }

            rms_a = sqrt(sum_a / audio_block_d.mags.size());
            rms_b = sqrt(sum_b / audio_block_c.mags.size());
            rms_c = sqrt(sum_c / audio_block_a.mags.size());
            rms_d = sqrt(sum_d / audio_block_b.mags.size());

            status = true;
        }
    }

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

    return status;
}

LiveDecoderSource* MorphGridModule::source() {
    return &my_source;
}
