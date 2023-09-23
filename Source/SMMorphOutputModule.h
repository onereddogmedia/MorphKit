// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_MORPH_OUTPUT_MODULE_HH
#define SPECTMORPH_MORPH_OUTPUT_MODULE_HH

#include "SMMorphOperatorModule.h"
#include "SMMorphOutput.h"
#include "SMMorphPlanVoice.h"

namespace SpectMorph {

class MorphOutputModule : public MorphOperatorModule {
    const struct MorphOutput::Config* cfg;
    std::vector<MorphOperatorModule*> out_ops;
    std::vector<class EffectDecoder*> out_decoders;

  public:
    MorphOutputModule(MorphPlanVoice* voice);
    ~MorphOutputModule();

    void set_config(const MorphOperatorConfig* op_cfg);
    void process(size_t n_samples, float** values, size_t n_ports, const float* freq_in = nullptr);
    void prepareToPlay(float mix_freq);
    void retrigger(int channel, float freq, int midi_velocity, bool onset);
    void release();
};

} // namespace SpectMorph

#endif
