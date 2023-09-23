// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_LIVEDECODER_HH
#define SPECTMORPH_LIVEDECODER_HH

#include "SMAlignedArray.h"
#include "SMLiveDecoderSource.h"
#include "SMPolyPhaseInter.h"
#include "SMSineDecoder.h"
#include "SMWavSet.h"
#include <functional>
#include <vector>

namespace SpectMorph {

class LiveDecoder {
    struct PartialState {
        float freq;
        float phase;
    };
    std::experimental::fixed_capacity_vector<PartialState, AudioBlock::block_size> pstate[2], *last_pstate;

    struct PortamentoState {
        std::experimental::fixed_capacity_vector<float, 32768> buffer;
        double pos;
        bool active;

        enum { DELTA = 32 };
    } portamento_state;

    Audio* audio;

    IFFTSynth* ifft_synth;
    LiveDecoderSource* source;
    PolyPhaseInter* pp_inter;

    double frame_step;
    size_t zero_values_at_start_scaled;
    size_t loop_start_scaled;
    size_t loop_end_scaled;
    int loop_point;
    float current_freq;
    float current_mix_freq;

    size_t have_samples;
    size_t block_size;
    size_t pos;
    double env_pos;
    size_t frame_idx;
    double original_samples_norm_factor;

    AlignedArray<float, 16>* sse_samples;

    // timing related
    double start_env_pos = 0;

    Audio::LoopType get_loop_type();

    void process_internal(size_t n_values, float* audio_out, float portamento_stretch);

    void portamento_grow(double end_pos, float portamento_stretch);
    void portamento_shrink();

    void process_portamento(size_t n_values, const float* freq_in, float* audio_out);
    LiveDecoder();

  public:
    LiveDecoder(LiveDecoderSource* source);
    ~LiveDecoder();

    void prepareToPlay(float mix_freq);

    void retrigger(int channel, float freq, int midi_velocity, bool onset);
    void process(size_t n_values, const float* freq_in, float* audio_out);

    double current_pos() const;
    double fundamental_note() const;

    static size_t compute_loop_frame_index(size_t index, Audio* audio);
};

} // namespace SpectMorph
#endif
