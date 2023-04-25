// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_MORPH_OPERATOR_MODULE_HH
#define SPECTMORPH_MORPH_OPERATOR_MODULE_HH

#include "SMLiveDecoderSource.h"
#include "SMMorphOperator.h"

#include <string>

namespace SpectMorph {

class MorphPlanVoice;

class MorphOperatorModule {
  protected:
    MorphPlanVoice* morph_plan_voice;
    std::vector<MorphOperatorModule*> m_dependencies;
    int m_update_value_tag;
    MorphOperator::PtrID m_ptr_id;

    void clear_dependencies();
    void add_dependency(MorphOperatorModule* dep_mod);

  public:
    MorphOperatorModule(MorphPlanVoice* voice);
    virtual ~MorphOperatorModule();

    virtual void set_config(const MorphOperatorConfig* op_cfg) = 0;
    virtual void prepareToPlay(float mix_freq);
    virtual LiveDecoderSource* source();
    virtual float value();

    const std::vector<MorphOperatorModule*>& dependencies() const;
    int& update_value_tag();
    void set_ptr_id(MorphOperator::PtrID ptr_id);

    static MorphOperatorModule* create(const std::string& type, MorphPlanVoice* voice);
};

} // namespace SpectMorph

#endif
