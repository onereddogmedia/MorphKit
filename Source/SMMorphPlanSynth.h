// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_MORPH_PLAN_SYNTH_HH
#define SPECTMORPH_MORPH_PLAN_SYNTH_HH

#include "SMMorphOperator.h"
#include "SMMorphPlan.h"
#include <map>
#include <memory>

namespace SpectMorph {

class MorphPlanVoice;
class MorphModuleSharedState;

class MorphPlanSynth {
  protected:
    std::vector<MorphPlanVoice*> voices;
    std::vector<std::string> m_last_update_ids;
    std::string m_last_plan_id;
    std::vector<MorphOperatorConfigP> m_active_configs;

    float m_mix_freq;

    WavSetRepo* m_wav_set_repo;

  public:
    struct Update {
        struct Op {
            MorphOperator::PtrID ptr_id;
            std::string type;
            MorphOperatorConfig* config = nullptr;
        };
        bool cheap = false; // cheap update: same set of operators
        std::vector<Op> ops;
        std::vector<MorphOperatorConfigP> new_configs;
        std::vector<MorphOperatorConfigP> old_configs;
    };
    typedef std::shared_ptr<Update> UpdateP;

    MorphPlanSynth(size_t n_voices);
    ~MorphPlanSynth();

    void reload() {
        if (m_wav_set_repo) {
            m_wav_set_repo->reload();
        }
    }

    UpdateP prepare_update(MorphPlanPtr new_plan);
    void apply_update(UpdateP update);

    MorphPlanVoice* voice(size_t i) const;

    void set_mix_freq(float mix_freq);
    float mix_freq() const;
    bool have_output() const;
};

} // namespace SpectMorph

#endif
