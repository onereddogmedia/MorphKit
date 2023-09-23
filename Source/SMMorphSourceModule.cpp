// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMMorphSourceModule.h"
#include "SMLeakDebugger.h"
#include "SMMorphPlan.h"
#include "SMMorphPlanVoice.h"
#include "SMMorphSource.h"
#include "glib.h"

using namespace SpectMorph;

using std::max;
using std::string;
using std::vector;

static LeakDebugger leak_debugger("SpectMorph::MorphSourceModule");

static float freq_to_note(float freq) {
    return 69.0f + 12.0f * logf(freq / 440.0f) / logf(2.0f);
}

SimpleWavSetSource::SimpleWavSetSource() : wav_set(nullptr), active_audio(nullptr) {
}

SimpleWavSetSource::~SimpleWavSetSource() {
    // managed elsewhere
    wav_set = nullptr;
    active_audio = nullptr;
}

void SimpleWavSetSource::set_wav_set(WavSetRepo* wave_set_repo, const string& path) {
    if (wave_set_repo) {
        WavSet* new_wav_set = wave_set_repo->get(path);
        if (new_wav_set != wav_set) {
            wav_set = new_wav_set;
            active_audio = nullptr;
        }
    } else {
        wav_set = nullptr;
        active_audio = nullptr;
    }
}

void SimpleWavSetSource::prepareToPlay(float /*mix_freq*/) {
}

void SimpleWavSetSource::retrigger(int channel, float freq, int midi_velocity, bool) {
    Audio* best_audio = nullptr;
    float best_diff = 1e10;

    if (wav_set) {
        float note = freq_to_note(freq);
        // TODO: use of deallocated memory!
        for (vector<WavSetWave>::iterator wi = wav_set->waves.begin(); wi != wav_set->waves.end(); wi++) {
            Audio* audio = wi->audio;
            if (audio && wi->channel == channel && wi->velocity_range_min <= midi_velocity &&
                wi->velocity_range_max >= midi_velocity) {
                float audio_note = freq_to_note(audio->fundamental_freq);
                if (fabs(audio_note - note) < best_diff) {
                    best_diff = fabs(audio_note - note);
                    best_audio = audio;
                }
            }
        }
    }
    active_audio = best_audio;
}

Audio* SimpleWavSetSource::audio() {
    return active_audio;
}

AudioBlock* SimpleWavSetSource::audio_block(size_t index) {
    if (active_audio) {
        if (active_audio->end > 0) {
            if (index < (size_t)active_audio->end && index < active_audio->contents.size())
                return &active_audio->contents[index];
        } else {
            if (index < active_audio->contents.size())
                return &active_audio->contents[index];
        }
    }
    return nullptr;
}

MorphSourceModule::MorphSourceModule(MorphPlanVoice* voice) : MorphOperatorModule(voice) {
    leak_debugger.add(this);
}

MorphSourceModule::~MorphSourceModule() {
    leak_debugger.del(this);
}

LiveDecoderSource* MorphSourceModule::source() {
    return &my_source;
}

void MorphSourceModule::set_config(const MorphOperatorConfig* op_cfg) {
    auto cfg = dynamic_cast<const MorphSource::Config*>(op_cfg);

    my_source.set_wav_set(cfg->wav_set_repo, cfg->path);
}
