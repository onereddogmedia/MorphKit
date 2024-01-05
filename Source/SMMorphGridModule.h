// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_MORPH_GRID_MODULE_HH
#define SPECTMORPH_MORPH_GRID_MODULE_HH

#include "SMLiveDecoder.h"
#include "SMMorphGrid.h"
#include "SMMorphOperatorModule.h"
#include "SMMorphSourceModule.h"
#include "SMWavSet.h"

#include <array>

namespace SpectMorph {

class MorphGridModule : public MorphOperatorModule {
  public:
    struct InputNode {
        MorphOperatorModule* mod;
        double delta_db;
        bool has_source;
        SimpleWavSetSource source;
    };

  private:
    const MorphGrid::Config* cfg = nullptr;

    struct InputNodeMatrix {
        static constexpr int MAX_DIM = 7;
        std::array<InputNode, MAX_DIM * MAX_DIM> data;

      public:
        InputNode& operator()(int x, int y) {
            assert(x < MAX_DIM && y < MAX_DIM);
            return data[(size_t)(x + y * MAX_DIM)];
        }
    } input_nodes;

    // output
    Audio audio;

    struct MySource : public LiveDecoderSource {
        MorphGridModule* module;

        void prepareToPlay(float mix_freq) override;
        void retrigger(int channel, float freq, int midi_velocity, bool onset) override;
        Audio* audio() override;
        bool rt_audio_block(size_t index, RTAudioBlock& out_block) override;
    } my_source;

  public:
    MorphGridModule(MorphPlanVoice* voice);
    ~MorphGridModule();

    void set_config(const MorphOperatorConfig* cfg);
    LiveDecoderSource* source();
};

} // namespace SpectMorph

#endif
