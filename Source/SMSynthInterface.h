// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_SYNTH_INTERFACE_HH
#define SPECTMORPH_SYNTH_INTERFACE_HH

#include "SMMidiSynth.h"
#include "SMMorphPlanSynth.h"
#include "SMProject.h"

namespace SpectMorph {

class SynthInterface {
    Project* m_project = nullptr;

  public:
    SynthInterface(Project* project) {
        m_project = project;
    }
    Project* get_project() {
        return m_project;
    }
    template <class DATA> void send_control_event(const std::function<void(Project*)>& func, DATA* data = nullptr) {
        m_project->synth_take_control_event(new InstFunc(func, [data]() { delete data; }));
    }

    void emit_apply_update(MorphPlanSynth::UpdateP update) {
        /* ownership of update is transferred to the event */
        struct EventData {
            MorphPlanSynth::UpdateP update;
        }* event_data = new EventData;

        event_data->update = update;
        send_control_event([=](Project* project) { project->midi_synth()->apply_update(event_data->update); },
                           event_data);
    }
};

} // namespace SpectMorph

#endif
