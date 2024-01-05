// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMMidiSynth.h"
#include "SMMorphOutputModule.h"

#include <assert.h>

using namespace SpectMorph;

MidiSynth::MidiSynth(size_t n_voices) : morph_plan_synth(n_voices) {
    voices.clear();
    voices.resize(n_voices);

    for (size_t i = 0; i < n_voices; i++) {
        voices[i].mp_voice = morph_plan_synth.voice(i);
    }
}

MorphPlanSynth::UpdateP MidiSynth::prepare_update(MorphPlan& plan) {
    return morph_plan_synth.prepare_update(plan);
}

void MidiSynth::apply_update(MorphPlanSynth::UpdateP update) {
    morph_plan_synth.apply_update(update);
}

void MidiSynth::set_mix_freq(double mix_freq) {
    m_mix_freq = mix_freq;
    morph_plan_synth.set_mix_freq((float)mix_freq);
}

double MidiSynth::mix_freq() const {
    return m_mix_freq;
}
