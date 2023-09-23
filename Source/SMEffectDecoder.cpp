// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMEffectDecoder.h"
#include "SMMorphOutput.h"
#include "SMMorphOutputModule.h"
#include "SMMorphUtils.h"

using namespace SpectMorph;

namespace SpectMorph {

class EffectDecoderSource : public LiveDecoderSource {
    LiveDecoderSource* source;
    Audio m_audio;
    float m_skip;

  public:
    explicit EffectDecoderSource(LiveDecoderSource* source);

    void prepareToPlay(float mix_freq) override;
    void retrigger(int channel, float freq, int midi_velocity, bool onset) override;
    Audio* audio() override;
    AudioBlock* audio_block(size_t index) override;

    void set_skip(float m_skip);
};

} // namespace SpectMorph

void EffectDecoderSource::prepareToPlay(float mix_freq) {
    source->prepareToPlay(mix_freq);
}

void EffectDecoderSource::retrigger(int channel, float freq, int midi_velocity, bool onset) {
    source->retrigger(channel, freq, midi_velocity, onset);
}

Audio* EffectDecoderSource::audio() {
    return &m_audio;
}

AudioBlock* EffectDecoderSource::audio_block(size_t index) {
    const double time_ms = index + m_skip; // 1ms frame step

    return MorphUtils::get_normalized_block_ptr(source, time_ms);
}

void EffectDecoderSource::set_skip(float skip) {
    m_skip = skip;
}

EffectDecoderSource::EffectDecoderSource(LiveDecoderSource* source_) : source(source_), m_skip(0) {
    m_audio.fundamental_freq = 440;
    m_audio.mix_freq = 48000;
    m_audio.frame_size_ms = 1;
    m_audio.frame_step_ms = 1;
    m_audio.zeropad = 4;
    m_audio.loop_type = Audio::LOOP_NONE;
}

EffectDecoder::EffectDecoder(LiveDecoderSource* source)
    : original_source(source), skip_source(new EffectDecoderSource(source)) {

    chain_decoder.reset(new LiveDecoder(original_source));
    use_skip_source = false;
}

EffectDecoder::~EffectDecoder() {
}

void EffectDecoder::set_config(const MorphOutput::Config*) {
    if (use_skip_source) // use original source (no skip)
    {
        chain_decoder.reset(new LiveDecoder(original_source));
        use_skip_source = false;
    }
}

void EffectDecoder::prepareToPlay(float mix_freq) {
    chain_decoder->prepareToPlay(mix_freq);
}

void EffectDecoder::retrigger(int channel, float freq, int midi_velocity, bool onset) {
    g_assert(chain_decoder);

    chain_decoder->retrigger(channel, freq, midi_velocity, onset);
}

void EffectDecoder::process(size_t n_values, const float* freq_in, float* audio_out) {
    g_assert(chain_decoder);

    chain_decoder->process(n_values, freq_in, audio_out);
}

void EffectDecoder::release() {
}
