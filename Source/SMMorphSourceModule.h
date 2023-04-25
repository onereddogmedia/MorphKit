// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_MORPH_SOURCE_MODULE_HH
#define SPECTMORPH_MORPH_SOURCE_MODULE_HH

#include "SMMorphOperatorModule.h"
#include "SMWavSet.h"

namespace SpectMorph {

class SimpleWavSetSource : public LiveDecoderSource {
  private:
    WavSet* wav_set;
    Audio* active_audio;

  public:
    SimpleWavSetSource();
    SimpleWavSetSource(const SimpleWavSetSource& other) : wav_set(other.wav_set), active_audio(other.active_audio) {
    }
    SimpleWavSetSource operator=(SimpleWavSetSource&& other) {
        if (this != &other) {
            wav_set = other.wav_set;
            active_audio = other.active_audio;
        }
        return *this;
    }
    ~SimpleWavSetSource() override;

    void set_wav_set(WavSetRepo* wave_set_repo, const std::string& path);
    void prepareToPlay(float mix_freq) override;

    void retrigger(int channel, float freq, int midi_velocity) override;
    Audio* audio() override;
    AudioBlock* audio_block(size_t index) override;
};

class MorphSourceModule : public MorphOperatorModule {
  protected:
    SimpleWavSetSource my_source;

  public:
    MorphSourceModule(MorphPlanVoice* voice);
    ~MorphSourceModule();

    void set_config(const MorphOperatorConfig* op_cfg);
    LiveDecoderSource* source();
};
} // namespace SpectMorph

#endif
