// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_MORPH_UTILS_HH
#define SPECTMORPH_MORPH_UTILS_HH

#include <algorithm>
#include <stdint.h>
#include <vector>

#include "SMAudio.h"
#include "SMLiveDecoder.h"
#include "SMRTMemory.h"

namespace SpectMorph {

namespace MorphUtils {

struct FreqState {
    float freq_f;
    int used;
};

bool find_match(float freq, const FreqState* freq_state, size_t freq_state_size, size_t* index);
void init_freq_state(const std::vector<uint16_t>& fint, FreqState* freq_state);
void init_freq_state(const RTVector<uint16_t>& fint, FreqState* freq_state);

AudioBlock* get_normalized_block_ptr(LiveDecoderSource* source, double time_ms);
bool get_normalized_block(LiveDecoderSource* source, double time_ms, RTAudioBlock& out_audio_block);

} // namespace MorphUtils

} // namespace SpectMorph

#endif
