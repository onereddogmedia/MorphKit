// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMSineDecoder.h"
#include "SMAlignedArray.h"
#include "SMAudio.h"
#include "SMFft.h"
#include "SMIFftSynth.h"
#include "SMMath.h"
#include "SMUtils.h"
#include <assert.h>
#include <math.h>
#include <stdio.h>

using SpectMorph::AudioBlock;
using SpectMorph::SineDecoder;
using std::vector;

/**
 * \brief Constructor setting up the various decoding parameters
 *
 * @param mix_freq_    sample rate to be used for reconstruction
 * @param frame_size_  frame size (in samples)
 * @param frame_step_  frame step (in samples)
 * @param mode_        selects decoding algorithm to be used
 */
SineDecoder::SineDecoder(double fundamental_freq_, double mix_freq_, size_t frame_size_, size_t frame_step_, Mode mode_)
    : mix_freq(mix_freq_), fundamental_freq(fundamental_freq_), frame_size(frame_size_), frame_step(frame_step_),
      mode(mode_) {
    ifft_synth = nullptr;
}

SineDecoder::~SineDecoder() {
    if (ifft_synth) {
        delete ifft_synth;
        ifft_synth = nullptr;
    }
}

/**
 * \brief Function which decodes a part of the signal.
 *
 * This needs two adjecant frames as arguments.
 *
 * @param block        the current frame (the frame to be decoded)
 * @param next_block   the frame after the current frame
 * @param window       the reconstruction window used for MODE_PHASE_SYNC_OVERLAP
 */
void SineDecoder::process(const AudioBlock& block, const AudioBlock& next_block, const vector<double>& window,
                          vector<float>& decoded_sines) {
    /* phase synchronous reconstruction (no loops) */
    if (mode == MODE_PHASE_SYNC_OVERLAP) {
        AlignedArray<float, 16> aligned_decoded_sines(frame_size);
        for (size_t i = 0; i < block.freqs.size(); i++) {
            const double SA = double(frame_step) / double(frame_size) * 2.0;
            const double mag_epsilon = 1e-8;

            VectorSinParams params;
            params.mag = block.mags_f(i) * SA;
            if (params.mag > mag_epsilon) {
                params.mix_freq = mix_freq;
                params.freq = block.freqs_f(i) * fundamental_freq;
                params.phase = block.phases_f(i);
                params.mode = VectorSinParams::ADD;

                fast_vector_sinf(params, &aligned_decoded_sines[0], &aligned_decoded_sines[frame_size]);
            }
        }
        for (size_t t = 0; t < frame_size; t++)
            decoded_sines[t] = aligned_decoded_sines[t] * (float)window[t];
        return;
    } else if (mode == MODE_PHASE_SYNC_OVERLAP_IFFT) {
        const size_t block_size = frame_size;

        if (!ifft_synth)
            ifft_synth = new IFFTSynth(block_size, mix_freq, IFFTSynth::WIN_HANNING);

        ifft_synth->clear_partials();
        for (size_t i = 0; i < block.freqs.size(); i++) {
            const double SA = double(frame_step) / double(frame_size) * 2.0;
            const double mag_epsilon = 1e-8;

            const double mag = block.mags_f(i) * SA;
            const double freq = block.freqs_f(i) * fundamental_freq;
            if (mag > mag_epsilon)
                ifft_synth->render_partial(freq, mag, block.phases_f(i));
        }
        ifft_synth->get_samples(&decoded_sines[0]);
        return;
    }

    zero_float_block(decoded_sines.size(), &decoded_sines[0]);

    /* phase distorted reconstruction */
    vector<float> freqs(block.freqs.size());
    vector<float> nfreqs(next_block.freqs.size());

    for (size_t i = 0; i < freqs.size(); i++)
        freqs[i] = (float)(block.freqs_f(i) * fundamental_freq);

    for (size_t i = 0; i < nfreqs.size(); i++)
        nfreqs[i] = (float)(next_block.freqs_f(i) * fundamental_freq);

    size_t todo = freqs.size() + nfreqs.size();

    synth_fixed_phase = next_synth_fixed_phase;
    synth_fixed_phase.resize(freqs.size());
    next_synth_fixed_phase.resize(nfreqs.size());

    const double SIN_AMP = 1.0;
    const bool TRACKING_SYNTH = true;
    while (todo) {
        double best_delta = 1e10;
        int best_i = 0, best_j = 0; /* init to get rid of gcc warning */
        for (size_t i = 0; i < freqs.size(); i++) {
            for (size_t j = 0; j < nfreqs.size(); j++) {
                double delta = fabs(freqs[i] - nfreqs[j]) / freqs[i];
                if (freqs[i] >= 0 && nfreqs[j] >= 0 && delta < best_delta && delta < 0.1) {
                    best_delta = delta;
                    best_i = (int)i;
                    best_j = (int)j;
                }
            }
        }
        if (best_delta < 0.1) {
            double freq = freqs[(size_t)best_i];
            freqs[(size_t)best_i] = -1;
            double nfreq = nfreqs[(size_t)best_j];
            nfreqs[(size_t)best_j] = -1;
            double mag = block.mags_f((size_t)best_i);
            double nmag = block.mags_f((size_t)best_j);

            // fprintf (stderr, "%f | %f ==> %f | %f\n", freq, mag, nfreq, nmag);
            assert(fabs(nfreq - freq) / freq < 0.1);

            double phase_delta = 2 * M_PI * freq / mix_freq;
            double nphase_delta = 2 * M_PI * nfreq / mix_freq;
            double phase = synth_fixed_phase[(size_t)best_i];
            if (TRACKING_SYNTH) {
                for (size_t i = 0; i < frame_step; i++) {
                    double inter = double(i) / double(frame_step);

                    decoded_sines[i] += (float)(sin(phase) * ((1 - inter) * mag + inter * nmag) * SIN_AMP);
                    phase += (1 - inter) * phase_delta + inter * nphase_delta;
                    while (phase > 2 * M_PI)
                        phase -= 2 * M_PI;
                }
                next_synth_fixed_phase[(size_t)best_j] = phase;
            } else {
                for (size_t i = 0; i < frame_size; i++) {
                    decoded_sines[i] += (float)(sin(phase) * window[i] * mag * SIN_AMP);
                    phase += phase_delta;
                    while (phase > 2 * M_PI)
                        phase -= 2 * M_PI;
                    // nfreq phase required -> ramp
                    if (i == frame_step - 1)
                        next_synth_fixed_phase[(size_t)best_j] = phase;
                }
            }
            todo -= 2;
        } else {
            for (size_t from = 0; from < freqs.size(); from++) {
                if (freqs[from] > -1) {
                    double freq = freqs[from];
                    freqs[from] = -1;
                    double mag = block.mags_f(from);

                    // fprintf (stderr, "%f | %f   >>> \n", freq, mag);

                    double phase_delta = 2 * M_PI * freq / mix_freq;
                    double phase = synth_fixed_phase[from];
                    if (TRACKING_SYNTH) {
                        for (size_t i = 0; i < frame_step; i++) {
                            double inter = double(i) / double(frame_step);

                            decoded_sines[i] += (float)(sin(phase) * (1 - inter) * mag * SIN_AMP);
                            phase += phase_delta;
                            while (phase > 2 * M_PI)
                                phase -= 2 * M_PI;
                        }
                    } else {
                        for (size_t i = 0; i < frame_size; i++) {
                            decoded_sines[i] += (float)(sin(phase) * window[i] * mag * SIN_AMP);
                            phase += phase_delta;
                            while (phase > 2 * M_PI)
                                phase -= 2 * M_PI;
                        }
                    }
                    todo--;
                }
            }
            for (size_t to = 0; to < nfreqs.size(); to++) {
                if (nfreqs[to] > -1) {
                    double freq = nfreqs[to];
                    nfreqs[to] = -1;
                    double mag = next_block.mags_f(to);

                    // fprintf (stderr, "%f | %f   <<< \n", freq, mag);

                    double phase_delta = 2 * M_PI * freq / mix_freq;
                    double phase = 0;
                    if (TRACKING_SYNTH) {
                        for (size_t i = 0; i < frame_step; i++) {
                            double inter = double(i) / double(frame_step);

                            decoded_sines[i] += (float)(sin(phase) * inter * mag * SIN_AMP); /* XXX */
                            phase += phase_delta;
                            while (phase > 2 * M_PI)
                                phase -= 2 * M_PI;
                        }
                        next_synth_fixed_phase[to] = phase;
                    }
                    todo--;
                }
            }
        }
    }
}
