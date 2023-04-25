// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_PROJECT_HH
#define SPECTMORPH_PROJECT_HH

#include "SMMorphPlan.h"
#include "SMObject.h"
#include "SMWavSet.h"

#include <mutex>
#include <thread>

namespace SpectMorph {

class MidiSynth;
class SynthInterface;
class MorphWavSource;

class SynthControlEvent {
  public:
    virtual void run_rt(Project* project) = 0;
    virtual ~SynthControlEvent() {
    }
};

struct InstFunc : public SynthControlEvent {
    std::function<void(Project*)> func;
    std::function<void()> free_func;

  public:
    InstFunc(const std::function<void(Project*)>& func_, const std::function<void()>& free_func_)
        : func(func_), free_func(free_func_) {
    }
    ~InstFunc() {
        free_func();
    }
    void run_rt(Project* project) {
        func(project);
    }
};

class ControlEventVector {
    std::vector<std::unique_ptr<SynthControlEvent>> events;
    bool clear = false;

  public:
    void take(SynthControlEvent* ev);
    void run_rt(Project* project);
};

class Project : public SignalReceiver {
  public:
    enum class StorageModel { COPY, REFERENCE };

  private:
    std::vector<std::shared_ptr<WavSet>> wav_sets;

    std::unique_ptr<MidiSynth> m_midi_synth;
    double m_mix_freq = 0;
    RefPtr<MorphPlan> m_morph_plan;
    std::vector<unsigned char> m_last_plan_data;
    bool m_state_changed_notify = false;

    std::mutex m_synth_mutex;
    ControlEventVector m_control_events;   // protected by synth mutex
    std::vector<std::string> m_out_events; // protected by synth mutex
    bool m_state_changed = false;          // protected by synth mutex

    std::unique_ptr<SynthInterface> m_synth_interface;

    std::vector<MorphWavSource*> list_wav_sources();

    void on_plan_changed();
    void on_operator_added(MorphOperator* op);
    void on_operator_removed(MorphOperator* op);

  public:
    Project();

    void add_rebuild_result(int object_id, WavSet* wav_set);
    void clear_wav_sets();

    std::shared_ptr<WavSet> get_wav_set(int object_id);

    void synth_take_control_event(SynthControlEvent* event);

    std::mutex& synth_mutex() {
        /* the synthesis thread will typically not block on synth_mutex
         * instead, it tries locking it, and if that fails, continues
         *
         * the ui thread will block on the synth_mutex to enqueue events,
         * parameter changes (in form of a new morph plan, volume, ...)
         * and so on
         *
         * if the synthesis thread can obtain a lock, it will then be
         * able to process these events to update its internal state
         * and also send notifications back to the ui
         */
        return m_synth_mutex;
    }
    bool try_update_synth();
    void set_mix_freq(double mix_freq);
    void set_state_changed_notify(bool notify);
    void state_changed();
    bool voices_active();

    std::vector<std::string> notify_take_events();
    SynthInterface* synth_interface() const;
    MidiSynth* midi_synth() const;
    MorphPlanPtr morph_plan() const;

    Signal<double> signal_volume_changed;
};

} // namespace SpectMorph

#endif
