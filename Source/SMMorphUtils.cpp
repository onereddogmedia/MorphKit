// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMMorphUtils.h"
#include "SMMath.h"

#include <algorithm>

using std::min;
using std::vector;

namespace SpectMorph {

namespace MorphUtils {

static bool fs_cmp(const FreqState& fs1, const FreqState& fs2) {
    return fs1.freq_f < fs2.freq_f;
}

bool find_match(float freq, const FreqState* freq_state, size_t freq_state_size, size_t* index) {
    const float freq_start = freq - 0.5f;
    const float freq_end = freq + 0.5f;

    double min_diff = 1e20;
    size_t best_index = 0; // initialized to avoid compiler warning

    FreqState start_freq_state = {freq_start, 0};
    const FreqState* start_ptr = std::lower_bound(freq_state, freq_state + freq_state_size, start_freq_state, fs_cmp);
    size_t i = (size_t)(start_ptr - freq_state);

    while (i < freq_state_size && freq_state[i].freq_f < freq_end) {
        if (!freq_state[i].used) {
            double diff = fabs(freq - freq_state[i].freq_f);
            if (diff < min_diff) {
                best_index = i;
                min_diff = diff;
            }
        }
        i++;
    }
    if (min_diff < 0.5) {
        *index = best_index;
        return true;
    }
    return false;
}

void init_freq_state(const vector<uint16_t>& fint, FreqState* freq_state) {
    for (size_t i = 0; i < fint.size(); i++) {
        freq_state[i].freq_f = (float)sm_ifreq2freq(fint[i]);
        freq_state[i].used = 0;
    }
}

void init_freq_state(const RTVector<uint16_t>& fint, FreqState* freq_state) {
    for (size_t i = 0; i < fint.size(); i++) {
        freq_state[i].freq_f = (float)sm_ifreq2freq(fint[i]);
        freq_state[i].used = 0;
    }
}

bool get_normalized_block(LiveDecoderSource* source, double time_ms, RTAudioBlock& out_audio_block) {
    if (!source)
        return false;

    Audio* audio = source->audio();
    if (!audio)
        return false;

    if (audio->loop_type == Audio::LOOP_TIME_FORWARD) {
        const double loop_start_ms = audio->loop_start * 1000.0 / audio->mix_freq;
        const double loop_end_ms = audio->loop_end * 1000.0 / audio->mix_freq;

        if (loop_start_ms >= loop_end_ms) {
            /* loop_start_index usually should be less than loop_end_index, this is just
             * to handle corner cases and pathological cases
             */
            time_ms = min(time_ms, loop_start_ms);
        } else if (time_ms > loop_end_ms) {
            /* compute loop position: ensure that time_ms is in [loop_start_ms, loop_end_ms] */
            time_ms = fmod(time_ms - loop_start_ms, loop_end_ms - loop_start_ms) + loop_start_ms;
        }
    }

    size_t source_index = (size_t)(sm_round_positive(time_ms / audio->frame_step_ms) + audio->start);
    if (audio->loop_type == Audio::LOOP_FRAME_FORWARD || audio->loop_type == Audio::LOOP_FRAME_PING_PONG) {
        source_index = LiveDecoder::compute_loop_frame_index(source_index, audio);
    }

    return source->rt_audio_block(source_index, out_audio_block);
}

} // namespace MorphUtils

} // namespace SpectMorph
