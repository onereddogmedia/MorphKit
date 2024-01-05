// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_EFFECTDECODER_HH
#define SPECTMORPH_EFFECTDECODER_HH

#include "SMLiveDecoder.h"
#include "SMLiveDecoderSource.h"
#include "SMMorphOutput.h"

#include <memory>

namespace SpectMorph {

class EffectDecoderSource;
class MorphOutputModule;
class EffectDecoder {
    MorphOutputModule* output_module;

    LiveDecoder chain_decoder;

    float current_freq = 440;

  public:
    EffectDecoder(MorphOutputModule* output_module = nullptr, float mix_freq = 0);
    ~EffectDecoder();

    void set_config(const MorphOutput::Config* cfg, LiveDecoderSource* source, float mix_freq);

    void prepareToPlay(float mix_freq);
    void retrigger(int channel, float freq, int midi_velocity, bool onset);
    void process(RTMemoryArea& rt_memory_area, size_t n_values, const float* freq_in, float* audio_out);
    void release();
    bool done();

    double time_offset_ms() const;
};

} // namespace SpectMorph

#endif
