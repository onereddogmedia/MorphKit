#ifndef SPECTMORPH_MORPHKIT_HH
#define SPECTMORPH_MORPHKIT_HH

#include "SMMorphGrid.h"
#include "SMMorphOutput.h"
#include "SMMorphOutputModule.h"
#include "SMMorphPlanSynth.h"
#include "SMMorphPlanVoice.h"
#include "SMProject.h"
#include "SMSynthInterface.h"

namespace SpectMorph {

const size_t kMaxBlockSize = 4096;

class MorphKit {
  public:
    MorphKit(const std::string& path, const std::string& icloud);
    ~MorphKit();

    void setSampleRate(const double newRate);
    void setModel(const uint model, const std::string& name);
    void setLinear(const float value);

    MorphPlanPtr morph_plan() {
        return project->morph_plan();
    }

    void wavSetReload();

  private:
    void makeMorphPlan();

    std::unique_ptr<SpectMorph::Project> project;
};

class MorphKitVoice {
  public:
    MorphKitVoice(MorphPlanPtr plan) : mp_voice(nullptr), morph_plan(plan) {
    }

    void startNote(const float note_frequency, const int8_t midi_velocity);
    void setControl(const double x, const double y, const double gain_a, const double gain_b, const double gain_c,
                    const double gain_d);
    void prepareToPlay(const size_t voice_id);

    template <typename sampleType> void render(const double frequency, sampleType* out, const size_t numSamples);

    void getRMS(float& rms_a, float& rms_b, float& rms_c, float& rms_d);

  private:
    MorphPlanVoice* mp_voice;
    MorphPlanPtr morph_plan;
};

class MorphKitEncoder {
  public:
    void encodeToFile(std::vector<float>& samples, int channels, float sampleRate, int bitsPerSample, int midi_note,
                      float fundamental_freq, int start, int end, Audio::LoopType loop_type, int loop_start,
                      int loop_end, bool auto_vol, float& norm_db, bool auto_pitch, const std::string& sample_name,
                      const std::string& path, std::function<bool()> kill_function);

  private:
    void applyAutoVolume(SpectMorph::Audio& audio);
    void applyAutoTune(SpectMorph::Audio& audio);
};

} // namespace SpectMorph

#endif
