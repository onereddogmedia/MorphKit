// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMIFftSynth.h"
#include "SMBlockUtils.h"
#include "SMMath.h"
#include <assert.h>
#include <stdio.h>

#include <map>
#include <mutex>

using namespace SpectMorph;

using std::map;
using std::vector;

static std::mutex table_mutex;
static map<size_t, IFFTSynthTable*> table_for_block_size;

namespace SpectMorph {
vector<float> IFFTSynth::sin_table;
}

IFFTSynth::IFFTSynth(size_t block_size_, double mix_freq_, WindowType win_type)
    : block_size(block_size_), mix_freq(mix_freq_) {
    std::lock_guard lg(table_mutex);

    zero_padding = 256;

    fft_work = nullptr;
    fft_workN = 0;

    table = table_for_block_size[block_size];
    if (!table) {
        const int range = 4;

        table = new IFFTSynthTable();

        const size_t win_size = block_size * (size_t)zero_padding;
        float* win = FFT::new_array_float(win_size);
        float* wspectrum = FFT::new_array_float(win_size);

        std::fill(win, win + win_size, 0); // most of it should be zero due to zeropadding
        for (size_t i = 0; i < block_size; i++) {
            if (i < block_size / 2)
                win[i] = (float)window_blackman_harris_92(double(block_size / 2 - i) / block_size * 2 - 1.0);
            else
                win[win_size - block_size + i] =
                    (float)window_blackman_harris_92(double(i - block_size / 2) / block_size * 2 - 1.0);
        }

        fft_work =
            FFT::fftar_float(fftar_float_plan, block_size * (size_t)zero_padding, win, wspectrum, fft_work, fft_workN);

        // compute complete (symmetric) expanded window transform for all frequency fractions
        for (int freq_frac = 0; freq_frac < zero_padding; freq_frac++) {
            for (int i = -range; i <= range; i++) {
                int pos = i * 256 - freq_frac;
                table->win_trans.push_back(wspectrum[abs(pos * 2)]);
            }
        }
        FFT::free_array_float(win);
        FFT::free_array_float(wspectrum);

        table->win_scale = FFT::new_array_float(block_size); // SSE
        for (size_t i = 0; i < block_size; i++)
            table->win_scale[(i + block_size / 2) % block_size] =
                (float)(window_cos(2.0 * i / block_size - 1.0) / window_blackman_harris_92(2.0 * i / block_size - 1.0));

        // we only need to do this once per block size (FIXME: not thread safe yet)
        table_for_block_size[block_size] = table;
    }
    if (sin_table.empty()) {
        // sin() table
        sin_table.resize(SIN_TABLE_SIZE);
        for (size_t i = 0; i < SIN_TABLE_SIZE; i++)
            sin_table[i] = (float)sin(i * 2 * M_PI / SIN_TABLE_SIZE);
    }

    if (win_type == WIN_BLACKMAN_HARRIS_92)
        win_scale = nullptr;
    else
        win_scale = table->win_scale;

    fft_in = FFT::new_array_float(block_size);
    fft_out = FFT::new_array_float(block_size);
    FFT::fftsr_destructive_float(fftsr_destructive_float_plan, block_size, fft_in, fft_out, nullptr);
    Block::zero((uint)block_size, fft_in);
    Block::zero((uint)block_size, fft_out);
    freq256_factor = 1 / mix_freq * block_size * zero_padding;
    mag_norm = 0.5 / block_size;
}

IFFTSynth::~IFFTSynth() {
    FFT::free_array_float(fft_in);
    FFT::free_array_float(fft_out);
    FFT::cleanup(fftar_float_plan);
    FFT::cleanup(fftsr_destructive_float_plan);
    FFT::free_array_float(fft_work);
    fft_work = nullptr;
    fft_workN = 0;
}

void IFFTSynth::get_samples(float* samples, OutputMode output_mode) {
    FFT::fftsr_destructive_float(fftsr_destructive_float_plan, block_size, fft_in, fft_out, fft_work);

    if (win_scale)
        Block::mul((uint)block_size, fft_out, win_scale);

    if (output_mode == REPLACE) {
        memcpy(samples, &fft_out[block_size / 2], sizeof(float) * block_size / 2);
        memcpy(&samples[block_size / 2], fft_out, sizeof(float) * block_size / 2);
    } else if (output_mode == ADD) {
        Block::add((uint)block_size / 2, samples, fft_out + block_size / 2);
        Block::add((uint)block_size / 2, samples + block_size / 2, fft_out);
    } else {
        assert(false);
    }
}

double IFFTSynth::quantized_freq(double mf_freq) {
    const int freq256 = sm_round_positive(mf_freq * freq256_factor);
    const double qfreq = freq256 * (1 / 256.0);
    const double mf_qfreq = qfreq / block_size * mix_freq;

    return mf_qfreq;
}

void IFFTSynth::precompute_tables() {
    // trigger fftw planning which can be slow
    FFT::fftsr_destructive_float(fftsr_destructive_float_plan, block_size, fft_in, fft_out, fft_work);
}
