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
    static constexpr size_t PARTIAL_STATE_RESERVE = 2048; // maximum number of partials to expect
    static constexpr size_t MAX_N_VALUES = 64;            // maximum number of values to process at once

    struct PartialState {
        float freq;
        float phase;
    };
    std::vector<PartialState> pstate[2], *last_pstate;

    struct PortamentoState {
        std::vector<float> buffer;
        double pos;
        bool active;

        enum { DELTA = 32 };
    } portamento_state;

    WavSet* smset;
    Audio* audio;

    size_t block_size;
    IFFTSynth* ifft_synth;
    LiveDecoderSource* source;
    PolyPhaseInter* pp_inter;
    RTMemoryArea* rt_memory_area = nullptr;

    double frame_step;
    size_t zero_values_at_start_scaled;
    size_t loop_start_scaled;
    size_t loop_end_scaled;
    int loop_point;
    float current_freq;
    float mix_freq;

    size_t have_samples;
    size_t pos;
    double env_pos;
    size_t frame_idx;

    AlignedArray<float, 16>* sse_samples;

    // timing related
    double start_env_pos = 0;
    bool in_process = false;

    Audio::LoopType get_loop_type();

    void process_internal(size_t n_values, float* audio_out, float portamento_stretch);

    void portamento_grow(double end_pos, float portamento_stretch);
    void portamento_shrink();

    void process_portamento(size_t n_values, const float* freq_in, float* audio_out);
    LiveDecoder();

  public:
    LiveDecoder(float mix_freq);
    LiveDecoder(WavSet* smset, float mix_freq);
    LiveDecoder(LiveDecoderSource* source, float mix_freq);
    ~LiveDecoder();

    void set_source(LiveDecoderSource* source);

    void prepareToPlay(float mix_freq);

    void retrigger(int channel, float freq, int midi_velocity, bool onset);
    void process(RTMemoryArea& rt_memory_area, size_t n_values, const float* freq_in, float* audio_out);

    double current_pos() const;
    double fundamental_note() const;

    static size_t compute_loop_frame_index(size_t index, Audio* audio);

    double time_offset_ms() const;
};

} // namespace SpectMorph
#endif
