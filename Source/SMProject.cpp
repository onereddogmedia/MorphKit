// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMProject.h"
#include "SMLiveDecoder.h"
#include "SMMemOut.h"
#include "SMMidiSynth.h"
#include "SMMorphOutput.h"
#include "SMMorphOutputModule.h"
#include "SMProject.h"
#include "SMSynthInterface.h"

using namespace SpectMorph;

using std::map;
using std::set;
using std::string;
using std::vector;

void ControlEventVector::take(SynthControlEvent* ev) {
    // we'd rather run destructors in non-rt part of the code
    if (clear) {
        events.clear();
        clear = false;
    }

    events.emplace_back(ev);
}

void ControlEventVector::run_rt(Project* project) {
    if (!clear) {
        for (const auto& ev : events)
            ev->run_rt(project);

        clear = true;
    }
}

void ControlEventVector::destroy_all_events() {
    events.clear();
    clear = false;
}

bool Project::try_update_synth() {
    bool state_changed = false;
    // handle synth updates (if locking is possible without blocking)
    //  - apply new parameters
    //  - process events
    if (m_synth_mutex.try_lock()) {
        m_control_events.run_rt(this);
        state_changed = m_state_changed;
        m_state_changed = false;

        m_synth_mutex.unlock();
    }
    return state_changed;
}

void Project::synth_take_control_event(SynthControlEvent* event) {
    std::lock_guard<std::mutex> lg(m_synth_mutex);
    m_control_events.take(event);
}

void Project::add_rebuild_result(int object_id, WavSet* wav_set) {
    // this function runs in audio thread
    int s = object_id + 1;
    if (s > (int)wav_sets.size())
        wav_sets.resize((size_t)s);

    wav_sets[(size_t)object_id] = std::shared_ptr<WavSet>(wav_set);
}

void Project::clear_wav_sets() {
    wav_sets.clear();
}

std::shared_ptr<WavSet> Project::get_wav_set(int object_id) {
    if (size_t(object_id) < wav_sets.size())
        return wav_sets[(size_t)object_id];
    else
        return nullptr;
}

SynthInterface* Project::synth_interface() const {
    return m_synth_interface.get();
}

MidiSynth* Project::midi_synth() const {
    return m_midi_synth.get();
}

Project::Project() : m_morph_plan(*this) {
    connect(m_morph_plan.signal_plan_changed, this, &Project::on_plan_changed);
    connect(m_morph_plan.signal_operator_added, this, &Project::on_operator_added);
    connect(m_morph_plan.signal_operator_removed, this, &Project::on_operator_removed);

    m_synth_interface.reset(new SynthInterface(this));

    /* avoid malloc in audio threads if wav sets are added */
    wav_sets.reserve(Project::WAV_SETS_RESERVE);
    m_midi_synth.reset(new MidiSynth(16));
}

void Project::set_mix_freq(double mix_freq) {
    /* if there are old control events, these cannot be executed anymore because we're
     * deleting the MidiSynth they refer to; we don't need a lock here because
     * this function is not supposed to be running while the audio thread is active
     */
    m_control_events.destroy_all_events();

    // not rt safe, needs to be called when synthesis thread is not running
    m_midi_synth->set_mix_freq(mix_freq);
    m_mix_freq = mix_freq;

    auto update = m_midi_synth->prepare_update(m_morph_plan);
    m_midi_synth->apply_update(update);
}

void Project::set_state_changed_notify(bool notify) {
    m_state_changed_notify = notify;
}

void Project::state_changed() {
    if (m_state_changed_notify)
        m_state_changed = true;
}

void Project::on_plan_changed() {
    state_changed();

    MorphPlanSynth::UpdateP update = m_midi_synth->prepare_update(m_morph_plan);
    m_synth_interface->emit_apply_update(update);
}

void Project::on_operator_added(MorphOperator*) {
}

void Project::on_operator_removed(MorphOperator*) {
}

MorphPlan* Project::morph_plan() {
    return &m_morph_plan;
}
