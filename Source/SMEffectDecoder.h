// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_EFFECTDECODER_HH
#define SPECTMORPH_EFFECTDECODER_HH

#include "SMLiveDecoder.h"
#include "SMLiveDecoderSource.h"
#include "SMMorphOutput.h"
#include "SMMorphOutputModule.h"

#include <memory>

namespace SpectMorph {

class EffectDecoderSource;
class EffectDecoder {
    LiveDecoderSource* original_source;

    bool use_skip_source;

    std::unique_ptr<EffectDecoderSource> skip_source;
    std::unique_ptr<LiveDecoder> chain_decoder;

  public:
    EffectDecoder(LiveDecoderSource* source);
    ~EffectDecoder();

    void set_config(const MorphOutput::Config* cfg);
    void prepareToPlay(float mix_freq);

    void retrigger(int channel, float freq, int midi_velocity, bool onset);
    void process(size_t n_values, const float* freq_in, float* audio_out);
    void release();
};

} // namespace SpectMorph

#endif
