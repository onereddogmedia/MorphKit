// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMWavData.h"
#include "SMMath.h"

// #include <sndfile.h>
#include <assert.h>

using namespace SpectMorph;

using std::string;
using std::vector;

WavData::WavData() {
    clear();
}

void WavData::clear() {
    m_samples.clear();

    m_n_channels = 0;
    m_bit_depth = 0;
    m_mix_freq = 0;
    m_error_blurb = "";
}

typedef int sf_count_t;

namespace {
struct VirtualData {
    vector<unsigned char>* mem = nullptr;
    sf_count_t offset = 0;
};
} // namespace

WavData::WavData(const vector<float>& samples, int n_channels, float mix_freq, int bit_depth) {
    m_samples = samples;
    m_n_channels = n_channels;
    m_mix_freq = mix_freq;
    m_bit_depth = bit_depth;
}

template <typename T> std::vector<T> slice(std::vector<T> const& v, int m, int n) {
    auto first = v.cbegin() + m;
    auto last = v.cbegin() + n + 1;

    std::vector<T> vec(first, last);
    return vec;
}

void WavData::crop(int start, int end) {
    m_samples = slice(m_samples, start, end);
}

void WavData::prepend(const vector<float>& samples) {
    assert(samples.size() % (size_t)m_n_channels == 0);

    m_samples.insert(m_samples.begin(), samples.begin(), samples.end());
}

float WavData::operator[](size_t pos) const {
    assert(pos < m_samples.size());

    return m_samples[pos];
}

float WavData::mix_freq() const {
    return m_mix_freq;
}

int WavData::n_channels() const {
    return m_n_channels;
}

int WavData::bit_depth() const {
    return m_bit_depth;
}

const vector<float>& WavData::samples() const {
    return m_samples;
}

size_t WavData::n_values() const {
    return m_samples.size();
}

const char* WavData::error_blurb() const {
    return m_error_blurb.c_str();
}
