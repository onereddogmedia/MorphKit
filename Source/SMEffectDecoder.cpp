// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMEffectDecoder.h"
#include "SMMorphOutput.h"
#include "SMMorphOutputModule.h"
#include "SMMorphUtils.h"

using namespace SpectMorph;

namespace SpectMorph {

class EffectDecoderSource : public LiveDecoderSource {
    LiveDecoderSource* m_source = nullptr;
    Audio m_audio;

  public:
    EffectDecoderSource();

    void prepareToPlay(float mix_freq) override;
    void retrigger(int channel, float freq, int midi_velocity, bool onset) override;
    Audio* audio() override;
    bool rt_audio_block(size_t index, RTAudioBlock& out_block) override;

    void set_source(LiveDecoderSource* source);
};

} // namespace SpectMorph

void EffectDecoderSource::prepareToPlay(float mix_freq) {
    if (m_source)
        m_source->prepareToPlay(mix_freq);
}

void EffectDecoderSource::retrigger(int channel, float freq, int midi_velocity, bool onset) {
    if (m_source)
        m_source->retrigger(channel, freq, midi_velocity, onset);
}

Audio* EffectDecoderSource::audio() {
    return &m_audio;
}

bool EffectDecoderSource::rt_audio_block(size_t index, RTAudioBlock& out_block) {
    const double time_ms = index; // 1ms frame step

    return MorphUtils::get_normalized_block(m_source, time_ms, out_block);
}

void EffectDecoderSource::set_source(LiveDecoderSource* source) {
    m_source = source;
}

EffectDecoderSource::EffectDecoderSource() {
    m_audio.fundamental_freq = 440;
    m_audio.mix_freq = 48000;
    m_audio.frame_size_ms = 1;
    m_audio.frame_step_ms = 1;
    m_audio.zeropad = 4;
    m_audio.loop_type = Audio::LOOP_NONE;
}

EffectDecoder::EffectDecoder(MorphOutputModule* output_module_, float mix_freq)
    : output_module(output_module_), chain_decoder(mix_freq) {
}

EffectDecoder::~EffectDecoder() {
}

void EffectDecoder::set_config(const MorphOutput::Config* /*cfg*/, LiveDecoderSource* source, float /*mix_freq*/) {
    chain_decoder.set_source(source);
}

void EffectDecoder::prepareToPlay(float mix_freq) {
    chain_decoder.prepareToPlay(mix_freq);
}

void EffectDecoder::retrigger(int channel, float freq, int midi_velocity, bool onset) {
    current_freq = freq;

    chain_decoder.retrigger(channel, freq, midi_velocity, onset);
}

void EffectDecoder::process(RTMemoryArea& rt_memory_area, size_t n_values, const float* freq_in, float* audio_out) {
    chain_decoder.process(rt_memory_area, n_values, freq_in, audio_out);
}

void EffectDecoder::release() {
}

bool EffectDecoder::done() {
    return false;
}

double EffectDecoder::time_offset_ms() const {
    return chain_decoder.time_offset_ms();
}
