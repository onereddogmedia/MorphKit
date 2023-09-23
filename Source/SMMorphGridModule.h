// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_MORPH_GRID_MODULE_HH
#define SPECTMORPH_MORPH_GRID_MODULE_HH

#include "SMMorphGrid.h"
#include "SMMorphOperatorModule.h"
#include "SMMorphSourceModule.h"
#include "SMWavSet.h"

namespace SpectMorph {

class MorphGridModule : public MorphOperatorModule {
  public:
    struct InputNode {
        bool has_source;
        SimpleWavSetSource source;
    };

  private:
    const MorphGrid::Config* cfg = nullptr;

    std::vector<std::vector<InputNode>> input_node;

    // output
    Audio audio;
    AudioBlock audio_block;

    struct MySource : public LiveDecoderSource {
        // temporary blocks for morphing:
        AudioBlock audio_block_a;
        AudioBlock audio_block_b;
        AudioBlock audio_block_c;
        AudioBlock audio_block_d;
        AudioBlock audio_block_ab;
        AudioBlock audio_block_cd;

        MorphGridModule* module;

        void prepareToPlay(float mix_freq) override;
        void retrigger(int channel, float freq, int midi_velocity, bool onset) override;
        Audio* audio() override;
        AudioBlock* audio_block(size_t index) override;
    } my_source;

  public:
    MorphGridModule(MorphPlanVoice* voice);
    ~MorphGridModule();

    void set_config(const MorphOperatorConfig* cfg);
    LiveDecoderSource* source();
};

} // namespace SpectMorph

#endif
