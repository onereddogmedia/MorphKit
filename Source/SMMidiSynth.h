// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_MIDI_SYNTH_HH
#define SPECTMORPH_MIDI_SYNTH_HH

#include "SMMorphPlanSynth.h"

namespace SpectMorph {

class MidiSynth {
    class Voice {
      public:
        MorphPlanVoice* mp_voice;

        Voice() : mp_voice(nullptr) {
        }
        Voice(const Voice& other) : mp_voice(other.mp_voice) {
        }
        Voice& operator=(Voice&& other) {
            if (this != &other) {
                mp_voice = other.mp_voice;
            }
            return *this;
        }
        ~Voice() {
            mp_voice = nullptr;
        }
    };

    MorphPlanSynth morph_plan_synth;

    std::vector<Voice> voices;
    double m_mix_freq;

  public:
    MidiSynth(size_t n_voices);

    MorphPlanSynth::UpdateP prepare_update(MorphPlanPtr plan);
    void apply_update(MorphPlanSynth::UpdateP update);

    void set_mix_freq(double mix_freq);
    double mix_freq() const;

    MorphPlanSynth* plan_synth() {
        return &morph_plan_synth;
    }
};

} // namespace SpectMorph

#endif /* SPECTMORPH_MIDI_SYNTH_HH */
