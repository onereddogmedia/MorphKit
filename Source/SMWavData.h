// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_WAVE_DATA_HH
#define SPECTMORPH_WAVE_DATA_HH

#include "SMUtils.h"

#include <functional>
#include <vector>

// #include <sndfile.h>

typedef void SNDFILE;
typedef void SF_INFO;

namespace SpectMorph {

class WavData {
  private:
    std::vector<float> m_samples;
    float m_mix_freq;
    int m_n_channels;
    int m_bit_depth;
    std::string m_error_blurb;

  public:
    WavData();
    WavData(const std::vector<float>& samples, int n_channels, float mix_freq, int bit_depth);

    void clear();
    void crop(int start, int end);
    void prepend(const std::vector<float>& samples);

    float mix_freq() const;
    int n_channels() const;
    size_t n_values() const;
    int bit_depth() const;
    const std::vector<float>& samples() const;
    const char* error_blurb() const;

    float operator[](size_t pos) const;
};

} // namespace SpectMorph

#endif /* SPECTMORPH_WAVE_DATA_HH */
