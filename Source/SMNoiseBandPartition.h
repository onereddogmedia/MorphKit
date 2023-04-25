// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_NOISE_BAND_PARTITION_HH
#define SPECTMORPH_NOISE_BAND_PARTITION_HH

#include <vector>

#include "SMAudio.h"
#include "SMRandom.h"
#include "glib.h"
#include <stdint.h>

namespace SpectMorph {

class NoiseBandPartition {
    std::vector<int> band_count;
    std::vector<int> band_start;

  public:
    NoiseBandPartition(double mix_freq);
    void noise_envelope_to_spectrum(SpectMorph::Random& random_gen, const AudioBlock::Block& envelope, float* spectrum,
                                    float scale);

    int bins_per_band(size_t band) {
        g_return_val_if_fail(band < band_count.size(), 0);

        return band_count[band];
    }

    constexpr static int n_bands = 32;
    constexpr static int n_spectrum_bins = 4096;
};

} // namespace SpectMorph

#endif
