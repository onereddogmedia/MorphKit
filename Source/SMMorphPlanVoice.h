// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_MORPH_PLAN_VOICE_HH
#define SPECTMORPH_MORPH_PLAN_VOICE_HH

#include "SMMorphOperatorModule.h"
#include "SMMorphPlan.h"
#include "SMMorphPlanSynth.h"

namespace SpectMorph {

class MorphOutputModule;
class MorphPlanSynth;

class MorphPlanVoice {
  protected:
    struct OpModule {
        MorphOperatorModule* module = nullptr;
        MorphOperator::PtrID ptr_id;
        MorphOperatorConfig* config = nullptr;
    };
    std::vector<OpModule> modules;

    std::vector<double> m_control_input;
    std::vector<double> m_control_output;
    MorphOutputModule* m_output;
    float m_mix_freq;
    MorphPlanSynth* m_morph_plan_synth;

    void clear_modules();
    void create_modules(MorphPlanSynth::UpdateP update);
    void configure_modules();

  public:
    MorphPlanVoice(MorphPlanSynth* synth);
    ~MorphPlanVoice();

    void cheap_update(MorphPlanSynth::UpdateP update);
    void full_update(MorphPlanSynth::UpdateP update);

    MorphOperatorModule* module(const MorphOperatorPtr& ptr);

    double control_input(double value, MorphOperator::ControlType ctype);
    void set_control_input(unsigned int i, double value);

    double control_output(unsigned int i);
    void set_control_output(unsigned int i, double value);

    void set_mix_freq(float freq);
    float mix_freq() const;

    MorphOutputModule* output();
    MorphPlanSynth* morph_plan_synth() const;
};

} // namespace SpectMorph

#endif
