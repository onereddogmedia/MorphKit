// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "SMMath.h"
#include "SMNoiseBandPartition.h"

using namespace SpectMorph;
using std::vector;

static double mel_to_hz(double mel) {
    return 700 * (exp(mel / 1127.0) - 1);
}

NoiseBandPartition::NoiseBandPartition(double mix_freq) : band_count((size_t)n_bands), band_start((size_t)n_bands) {
    int d = 0;
    /* assign each d to a band */
    vector<int> band_from_d((size_t)n_spectrum_bins);
    std::fill(band_from_d.begin(), band_from_d.end(), -1);
    for (int band = 0; band < n_bands; band++) {
        double mel_low = 30 + 4000.0 / n_bands * band;
        double mel_high = 30 + 4000.0 / n_bands * (band + 1);
        double hz_low = mel_to_hz(mel_low);
        double hz_high = mel_to_hz(mel_high);

        /* skip frequencies which are too low to be in lowest band */
        double f_hz = mix_freq / 2.0 * d / n_spectrum_bins;
        if (band == 0) {
            while (f_hz < hz_low) {
                d += 2;
                f_hz = mix_freq / 2.0 * d / n_spectrum_bins;
            }
        }
        while (f_hz < hz_high && d < n_spectrum_bins) {
            if (d < (int)band_from_d.size()) {
                band_from_d[(size_t)d] = band;
                band_from_d[(size_t)(d + 1)] = band;
            }
            d += 2;
            f_hz = mix_freq / 2.0 * d / n_spectrum_bins;
        }
    }
    /* count bins per band */
    for (d = 0; d < n_spectrum_bins; d += 2) {
        int b = band_from_d[(size_t)d];
        if (b != -1) {
            assert(b >= 0 && b < int(n_bands));
            if (band_count[(size_t)b] == 0)
                band_start[(size_t)b] = d;
            band_count[(size_t)b]++;
        }
    }
}

void NoiseBandPartition::noise_envelope_to_spectrum(Random& random_gen, const AudioBlock::Block& envelope,
                                                    float* spectrum, float scale) {
    assert(envelope.size() == n_bands);

    uint32 random_data[(n_spectrum_bins + 7) / 8];

    random_gen.random_block((size_t)(n_spectrum_bins + 7) / 8, random_data);

    zero_float_block((size_t)n_spectrum_bins, spectrum);

    const uint8* random_data_byte = reinterpret_cast<uint8*>(&random_data[0]);

    for (size_t b = 0; b < n_bands; b++) {
        const float value = sm_idb2factor(envelope[b]) * scale;

        int start = band_start[b];
        int end = start + band_count[b] * 2;
        for (int d = start; d < end; d += 2) {
            /* Generate complex number with:
             *  - phase:     r / 256.0 * 2 * M_PI
             *  - magnitude: value
             */
            const uint8 r = random_data_byte[d / 2];

            spectrum[d] = int_cosf(r) * value;
            spectrum[d + 1] = int_sinf(r) * value;
        }
    }
}
