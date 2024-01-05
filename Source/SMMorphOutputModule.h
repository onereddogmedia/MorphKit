// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_MORPH_OUTPUT_MODULE_HH
#define SPECTMORPH_MORPH_OUTPUT_MODULE_HH

#include "SMEffectDecoder.h"
#include "SMMorphOperatorModule.h"
#include "SMMorphPlanVoice.h"
#include "SMRTMemory.h"

namespace SpectMorph {

class MorphOutputModule : public MorphOperatorModule {
    const MorphOutput::Config* cfg = nullptr;
    const TimeInfoGenerator* time_info_gen = nullptr;
    RTMemoryArea* m_rt_memory_area = nullptr;
    EffectDecoder decoder;

  public:
    MorphOutputModule(MorphPlanVoice* voice);
    ~MorphOutputModule();

    void set_config(const MorphOperatorConfig* op_cfg);
    void process(const TimeInfoGenerator& time_info, RTMemoryArea& rt_memory_area, size_t n_samples, float** values,
                 const float* freq_in = nullptr);
    void prepareToPlay(float mix_freq);
    void retrigger(const TimeInfo& time_info, int channel, float freq, int midi_velocity, bool onset);
    void release();
    bool done();
    TimeInfo compute_time_info() const;
    RTMemoryArea* rt_memory_area() const;
};

} // namespace SpectMorph

#endif
