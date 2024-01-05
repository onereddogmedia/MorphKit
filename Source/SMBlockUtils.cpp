// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMBlockUtils.h"

using namespace SpectMorph;

void Block::add(uint n_values, float* ovalues, const float* ivalues) {
    // auto vectorization for simple cases is quite good these days,
    // so we don't provide an SSE intrinsics implementation

    for (uint i = 0; i < n_values; i++)
        ovalues[i] += ivalues[i];
}

void Block::mul(uint n_values, float* ovalues, const float* ivalues) {
    // auto vectorization for simple cases is quite good these days,
    // so we don't provide an SSE intrinsics implementation

    for (uint i = 0; i < n_values; i++)
        ovalues[i] *= ivalues[i];
}

void Block::zero(uint n_values, float* ovalues) {
    for (uint i = 0; i < n_values; i++)
        ovalues[i] = 0;
}

void Block::range(uint n_values, const float* ivalues, float& min_value, float& max_value) {
    float minv, maxv;
    if (n_values) {
        minv = maxv = ivalues[0];

        for (uint i = 1; i < n_values; i++) {
            if (ivalues[i] < minv)
                minv = ivalues[i];
            if (ivalues[i] > maxv)
                maxv = ivalues[i];
        }
    } else {
        minv = maxv = 0;
    }
    min_value = minv;
    max_value = maxv;
}
