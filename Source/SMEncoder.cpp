// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMEncoder.h"
#include "SMAlignedArray.h"
#include "SMBlockUtils.h"
#include "SMFft.h"
#include "SMMath.h"
#include "SMRandom.h"
#include "SMUtils.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include <algorithm>
#include <cinttypes>
#include <complex>
#include <map>
#include <memory>

using namespace SpectMorph;
using std::complex;
using std::map;
using std::max;
using std::string;
using std::vector;

#define MIN(a, b) b

static double magnitude(vector<float>::iterator i) {
    return sqrt(*i * *i + *(i + 1) * *(i + 1));
}

// wraps phase in range [0:2*pi]
static double normalize_phase(double phase) {
    const double inv_2pi = 1.0 / (2.0 * M_PI);
    phase *= inv_2pi;
    phase -= floor(phase);
    return phase * (2.0 * M_PI);
}

#define debug(...) SpectMorph::Debug::debug("encoder", __VA_ARGS__)

EncoderParams::EncoderParams()
    : param_name_d({"peak-width", "min-frame-periods", "min-frame-size"}), param_name_s({"window"}) {
}

bool EncoderParams::add_config_entry(const string& param, const string& value) {
    for (auto name : param_name_d) // double parameters
    {
        if (name == param) {
            param_value_d[name] = sm_atof(value.c_str());
            return true;
        }
    }
    for (auto name : param_name_s) // string parameters
    {
        if (name == param) {
            param_value_s[name] = value;
            return true;
        }
    }
    return false; // unsupported parameter
}

bool EncoderParams::get_param(const string& param, double& value) const {
    if (find(param_name_d.begin(), param_name_d.end(), param) == param_name_d.end()) {
        fprintf(stderr, "error: encoder parameter '%s' was not defined\n", param.c_str());
        return false;
    }

    map<string, double>::const_iterator pi = param_value_d.find(param);
    if (pi == param_value_d.end()) {
        return false; /* not defined */
    } else {
        value = pi->second;
        return true;
    }
}

bool EncoderParams::get_param(const string& param, string& value) const {
    if (find(param_name_s.begin(), param_name_s.end(), param) == param_name_s.end()) {
        fprintf(stderr, "error: encoder parameter '%s' was not defined\n", param.c_str());
        return false;
    }

    map<string, string>::const_iterator pi = param_value_s.find(param);
    if (pi == param_value_s.end()) {
        return false; /* not defined */
    } else {
        value = pi->second;
        return true;
    }
}

static size_t make_odd(size_t n) {
    if (n & 1)
        return n;
    return n - 1;
}

void EncoderParams::setup_params(const WavData& wav_data, double new_fundamental_freq) {
    mix_freq = wav_data.mix_freq();
    zeropad = 4;
    fundamental_freq = new_fundamental_freq;

    // --- frame size & step ---
    double min_frame_periods, min_frame_size;
    if (!get_param("min-frame-periods", min_frame_periods))
        min_frame_periods = 4; // default: at least 4 periods of the fundamental per frame
    if (!get_param("min-frame-size", min_frame_size))
        min_frame_size = 40; // default: at least 40ms frames

    frame_size_ms = (float)min_frame_size;
    frame_size_ms = max<float>(frame_size_ms, 1000.0f / (float)fundamental_freq * (float)min_frame_periods);
    frame_step_ms = frame_size_ms / 4.0f;

    // --- convert block sizes in ms to sample counts ----
    frame_size = make_odd((size_t)(mix_freq * 0.001f * frame_size_ms));
    frame_step = (size_t)(mix_freq * 0.001f * frame_step_ms);

    /* compute block size from frame size (smallest 2^k value >= frame_size) */
    block_size = 1;
    while (block_size < frame_size)
        block_size *= 2;

    /* compute encoder window */
    window.resize(block_size);

    for (size_t i = 0; i < window.size(); i++) {
        if (i < frame_size)
            window[i] = (float)window_cos(2.0 * i / (frame_size - 1) - 1.0);
        else
            window[i] = 0;
    }
}

void EncoderParams::set_kill_function(const std::function<bool()>& new_kill_function) {
    kill_function = new_kill_function;
}

/**
 * Constructor which initializes the Encoders parameters.
 */
Encoder::Encoder(const EncoderParams& enc_params_) {
    assert(enc_params_.mix_freq > 0);
    assert(enc_params_.frame_step_ms > 0);
    assert(enc_params_.frame_size_ms > 0);
    assert(enc_params_.zeropad > 0);
    assert(enc_params_.frame_step > 0);
    assert(enc_params_.frame_size > 0);
    assert(enc_params_.block_size > 0);
    assert(enc_params_.fundamental_freq > 0);
    assert(enc_params_.window.size() == enc_params_.block_size);

    this->enc_params = enc_params_;

    loop_start = -1;
    loop_end = -1;
    loop_type = Audio::LOOP_NONE;
}

/**
 * This function computes the short-time-fourier-transform (STFT) of the input
 * signal using a window to cut the individual frames out of the sample.
 */
void Encoder::compute_stft(const WavData& multi_channel_wav_data, int channel) {
    /* deinterleave multi channel signal */
    vector<float> single_channel_signal;

    const int n_channels = multi_channel_wav_data.n_channels();
    for (int i = channel; i < (int)multi_channel_wav_data.n_values(); i += n_channels)
        single_channel_signal.push_back(multi_channel_wav_data[(size_t)i]);

    original_samples = single_channel_signal;

    WavData wav_data(single_channel_signal, 1, multi_channel_wav_data.mix_freq(), multi_channel_wav_data.bit_depth());

    /* encode single channel */
    zero_values_at_start = enc_params.frame_size - enc_params.frame_step / 2;
    vector<float> zero_values(zero_values_at_start);

    wav_data.prepend(zero_values);

    const uint64 n_values = wav_data.n_values();
    const size_t frame_size = enc_params.frame_size;
    const size_t block_size = enc_params.block_size;
    const size_t zeropad = enc_params.zeropad;
    const auto& window = enc_params.window;

    sample_count = n_values;

    vector<double> in(block_size * zeropad), out(block_size * zeropad + 2);

    float* fft_in = FFT::new_array_float(in.size());
    float* fft_out = FFT::new_array_float(in.size());
    float* fft_work = nullptr;
    size_t fft_workN = 0;
    FFT::PlanMap fftar_float_plan;

    for (uint64 pos = 0; pos < n_values; pos += enc_params.frame_step) {
        EncoderBlock audio_block;

        /* start with zero block, so the incomplete blocks at end are zeropadded */
        vector<float> block(block_size);

        for (size_t offset = 0; offset < block.size(); offset++) {
            if (pos + offset < wav_data.n_values())
                block[offset] = wav_data[pos + offset];
        }
        vector<float> debug_samples(block.begin(), block.end());
        Block::mul((uint)enc_params.block_size, &block[0], &window[0]);

        size_t j = in.size() - enc_params.frame_size / 2;
        for (vector<float>::const_iterator i = block.begin(); i != block.end(); i++)
            in[((j++) % in.size())] = *i;

        std::copy(in.begin(), in.end(), fft_in);
        fft_work = FFT::fftar_float(fftar_float_plan, in.size(), fft_in, fft_out, fft_work, fft_workN);
        std::copy(fft_out, fft_out + in.size(), out.begin());

        out[block_size * zeropad] = out[1];
        out[block_size * zeropad + 1] = 0;
        out[1] = 0;

        audio_block.noise.assign(out.begin(), out.end()); // <- will be overwritten by noise spectrum later on
        audio_block.debug_samples.assign(debug_samples.begin(), debug_samples.begin() + (long)frame_size);
        audio_blocks.push_back(audio_block);

        if (killed("_stft", audio_blocks.size() & 63))
            break; // break to avoid leaking fft_in, fft_out
    }
    FFT::free_array_float(fft_in);
    FFT::free_array_float(fft_out);
    FFT::free_array_float(fft_work);
    FFT::cleanup(fftar_float_plan);
}

namespace {

class QInterpolator {
    double a, b, c;

  public:
    QInterpolator(double y1, double y2, double y3) {
        a = (y1 + y3 - 2 * y2) / 2;
        b = (y3 - y1) / 2;
        c = y2;
    }
    double eval(double x) {
        return a * x * x + b * x + c;
    }
    double x_max() {
        return -b / (2 * a);
    }
};

} // namespace

/**
 * This function searches for peaks in the frame ffts. These are stored in frame_tracksels.
 */
void Encoder::search_local_maxima() {
    const size_t block_size = enc_params.block_size;
    const size_t frame_size = enc_params.frame_size;
    const size_t zeropad = enc_params.zeropad;
    const double mix_freq = enc_params.mix_freq;
    const auto& window = enc_params.window;

    // figure out normalization for window
    double window_weight = 0;
    for (size_t i = 0; i < frame_size; i++)
        window_weight += window[i];
    const double window_scale = 2.0 / window_weight;

    // initialize tracksel structure
    frame_tracksels.clear();
    frame_tracksels.resize(audio_blocks.size());

    // find maximum of all values
    double max_mag = 0;
    for (size_t n = 0; n < audio_blocks.size(); n++) {
        for (size_t d = 2; d < block_size * zeropad; d += 2) {
            max_mag = max(max_mag, magnitude(audio_blocks[n].noise.begin() + (long)d));
        }
    }

    for (size_t n = 0; n < audio_blocks.size(); n++) {
        vector<double> mag_values(audio_blocks[n].noise.size() / 2);
        for (size_t d = 0; d < block_size * zeropad; d += 2)
            mag_values[d / 2] = magnitude(audio_blocks[n].noise.begin() + (long)d);

        for (size_t d = 2; d < block_size * zeropad; d += 2) {
#if 0
	  double phase = atan2 (*(audio_blocks[n]->noise.begin() + d),
	                        *(audio_blocks[n]->noise.begin() + d + 1)) / 2 / M_PI;  /* range [-0.5 .. 0.5] */
#endif
            enum { PEAK_NONE, PEAK_SINGLE, PEAK_DOUBLE } peak_type = PEAK_NONE;

            if (mag_values[d / 2] > mag_values[d / 2 - 1] &&
                mag_values[d / 2] > mag_values[d / 2 + 1]) /* search for peaks in fft magnitudes */
            {
                /* single peak is the common case, where the magnitude of the middle value is
                 * larger than the magnitude of the left and right neighbour
                 */
                peak_type = PEAK_SINGLE;
            } else {
                double epsilon_fact = 1.0 + 1e-8;
                if (mag_values[d / 2] < mag_values[d / 2 + 1] * epsilon_fact &&
                    mag_values[d / 2] * epsilon_fact > mag_values[d / 2 + 1] &&
                    mag_values[d / 2] > mag_values[d / 2 - 1] && mag_values[d / 2] > mag_values[d / 2 + 2]) {
                    /* double peak is a special case, where two values in the spectrum have (almost) equal magnitude
                     * in this case, this magnitude must be larger than the value left and right of the _two_
                     * maximal values in the spectrum
                     */
                    peak_type = PEAK_DOUBLE;
                }
            }

            const double mag2 = db_from_factor(mag_values[d / 2] / max_mag, -100);

            if (peak_type != PEAK_NONE) {
                if (mag2 > -90) {
                    size_t ds, de;
                    for (ds = d / 2 - 1; ds > 0 && mag_values[ds] < mag_values[ds + 1]; ds--)
                        ;
                    for (de = d / 2 + 1; de < (mag_values.size() - 1) && mag_values[de] > mag_values[de + 1]; de++)
                        ;

                    const double normalized_peak_width = (de - ds) * frame_size / double(block_size * zeropad);

                    bool peak_ok;
                    double value;
                    if (enc_params.get_param("peak-width", value))
                        peak_ok = normalized_peak_width > value;
                    else
                        peak_ok = normalized_peak_width > 2.9;

                    if (peak_ok) {
                        const double mag1 = db_from_factor(mag_values[d / 2 - 1] / max_mag, -100);
                        const double mag3 = db_from_factor(mag_values[d / 2 + 1] / max_mag, -100);
                        // double freq = d / 2 * mix_freq / (block_size * zeropad); /* bin frequency */

                        QInterpolator mag_interp(mag1, mag2, mag3);
                        double x_max = mag_interp.x_max();
                        double tfreq = (d / 2 + x_max) * mix_freq / (block_size * zeropad);

                        double peak_mag_db = mag_interp.eval(x_max);
                        double peak_mag = db_to_factor(peak_mag_db) * max_mag;

                        // use the interpolation formula for the complex values to find the phase
                        QInterpolator re_interp(audio_blocks[n].noise[d - 2], audio_blocks[n].noise[d],
                                                audio_blocks[n].noise[d + 2]);
                        QInterpolator im_interp(audio_blocks[n].noise[d - 1], audio_blocks[n].noise[d + 1],
                                                audio_blocks[n].noise[d + 3]);
                        /*
                                          if (mag2 > -20)
                                            printf ("%f %f %f %f %f\n", phase, last_phase[d], phase_diff, phase_diff *
                           mix_freq / (block_size * zeropad) * overlap, tfreq);
                        */
                        Tracksel tracksel;
                        tracksel.frame = n;
                        tracksel.d = d;
                        tracksel.freq = tfreq;
                        tracksel.mag = peak_mag * window_scale;
                        tracksel.mag2 = mag2;
                        tracksel.next = nullptr;
                        tracksel.prev = nullptr;

                        const double re_mag = re_interp.eval(x_max);
                        const double im_mag = im_interp.eval(x_max);
                        double phase = atan2(im_mag, re_mag) + 0.5 * M_PI;
                        // correct for the odd-centered analysis
                        {
                            phase -= (frame_size - 1) / 2.0 / mix_freq * tracksel.freq * 2 * M_PI;
                            phase = normalize_phase(phase);
                        }
                        tracksel.phase = phase;

                        // FIXME: need a different criterion here
                        // mag2 > -30 doesn't track all partials
                        // mag2 > -60 tracks lots of junk, too
                        if (mag2 > -90 && tracksel.freq > 10)
                            frame_tracksels[n].push_back(tracksel);

                        if (peak_type == PEAK_DOUBLE)
                            d += 2;
                    }
                }
            }
        }

        if (killed("_maxima", n & 15))
            return;
    }
}

/// @cond
struct PeakIndex {
    double freq;
    vector<Tracksel>::iterator i;
    PeakIndex* prev;
    double prev_delta;

    PeakIndex(double freq_, vector<Tracksel>::iterator i_) : freq(freq_), i(i_), prev(nullptr), prev_delta(0) {
    }
};
/// @endcond

static bool partial_index_cmp(const PeakIndex& a, const PeakIndex& b) {
    return a.freq < b.freq;
}

/**
 * This function links the spectral peaks (contained in the Tracksel structure)
 * of successive frames together by setting the prev and next pointers. It
 * tries to minimize the frequency difference between the peaks that are linked
 * together, while using a threshold of 5% frequency derivation.
 */
void Encoder::link_partials() {
    for (size_t n = 0; n + 1 < audio_blocks.size(); n++) {
        // build sorted index for this frame
        vector<PeakIndex> current_index;
        for (vector<Tracksel>::iterator i = frame_tracksels[n].begin(); i != frame_tracksels[n].end(); i++)
            current_index.push_back(PeakIndex(i->freq, i));
        sort(current_index.begin(), current_index.end(), partial_index_cmp);

        // build sorted index for next frame
        vector<PeakIndex> next_index;
        for (vector<Tracksel>::iterator i = frame_tracksels[n + 1].begin(); i != frame_tracksels[n + 1].end(); i++)
            next_index.push_back(PeakIndex(i->freq, i));
        sort(next_index.begin(), next_index.end(), partial_index_cmp);

        vector<PeakIndex>::iterator ci = current_index.begin();
        vector<PeakIndex>::iterator ni = next_index.begin();
        if (ni != next_index.end()) // if current or next frame are empty (no peaks) there is nothing to do
        {
            while (ci != current_index.end()) {
                /*
                 * increment ni as long as incrementing it makes ni point to a
                 * better (closer) peak below ci's frequency
                 */
                vector<PeakIndex>::iterator inc_ni;
                do {
                    inc_ni = ni + 1;
                    if (inc_ni < next_index.end() && inc_ni->freq < ci->freq)
                        ni = inc_ni;
                } while (ni == inc_ni);

                /*
                 * possible candidates for a match are
                 * - ni      - which contains the greatest peak with a smaller frequency than ci->freq
                 * - ni + 1  - which contains the smallest peak with a greater frequency that ci->freq
                 * => choose the candidate which is closer to ci->freq
                 */
                vector<PeakIndex>::iterator besti = ni;
                if (ni + 1 < next_index.end() && fabs(ci->freq - (ni + 1)->freq) < fabs(ci->freq - ni->freq))
                    besti = ni + 1;

                const double delta = fabs(ci->freq - besti->freq) / ci->freq;
                if (delta < 0.05) /* less than 5% frequency derivation */
                {
                    if (!besti->prev || besti->prev_delta > delta) {
                        besti->prev = &(*ci);
                        besti->prev_delta = delta;
                    }
                }
                ci++;
            }

            /* link best matches (with the smallest frequency derivation) */
            for (ni = next_index.begin(); ni != next_index.end(); ni++) {
                if (ni->prev) {
                    Tracksel* crosslink_a = &(*ni->prev->i);
                    Tracksel* crosslink_b = &(*ni->i);
                    crosslink_a->next = crosslink_b;
                    crosslink_b->prev = crosslink_a;
                }
            }
        }
    }
}

/**
 * This function validates that the partials found by the peak linking have
 * good quality.
 */
void Encoder::validate_partials() {
    map<Tracksel*, bool> processed_tracksel;
    for (size_t n = 0; n < audio_blocks.size(); n++) {
        vector<Tracksel>::iterator i;
        for (i = frame_tracksels[n].begin(); i != frame_tracksels[n].end(); i++) {
            if (!processed_tracksel[&(*i)]) {
                double biggest_mag = -100;
                for (Tracksel* t = &(*i); t; t = t->next) {
                    biggest_mag = max(biggest_mag, t->mag2);
                    processed_tracksel[t] = true;
                }
                if (biggest_mag > -90) {
                    for (Tracksel* t = &(*i); t; t = t->next) {
                        audio_blocks[t->frame].freqs.push_back((float)t->freq);
                        audio_blocks[t->frame].mags.push_back((float)t->mag);
                        audio_blocks[t->frame].phases.push_back((float)t->phase);
                    }
                }
            }
        }
        if (killed("_validate", n & 63))
            return;
    }
}

/**
 * This function subtracts the partials from the audio signal, to get the
 * residue (remaining energy not corresponding to sine frequencies).
 */
void Encoder::spectral_subtract() {
    const size_t block_size = enc_params.block_size;
    const size_t frame_size = enc_params.frame_size;
    const size_t zeropad = (size_t)enc_params.zeropad;
    const auto& window = enc_params.window;

    float* fft_in = FFT::new_array_float(block_size * zeropad);
    float* fft_out = FFT::new_array_float(block_size * zeropad);
    float* fft_work = nullptr;
    size_t fft_workN = 0;
    FFT::PlanMap fftar_float_plan;

    for (uint64 frame = 0; frame < audio_blocks.size(); frame++) {
        AlignedArray<float, 16> signal(frame_size);
        for (size_t i = 0; i < audio_blocks[frame].freqs.size(); i++) {
            const double freq = audio_blocks[frame].freqs[i];
            const double mag = audio_blocks[frame].mags[i];
            const double phase = audio_blocks[frame].phases[i];

            VectorSinParams params;
            params.mix_freq = enc_params.mix_freq;
            params.freq = freq;
            params.phase = phase;
            params.mag = mag;
            params.mode = VectorSinParams::ADD;

            fast_vector_sinf(params, &signal[0], &signal[frame_size]);
        }
        vector<double> out(block_size * zeropad + 2);
        // apply window
        std::fill(fft_in, fft_in + block_size * zeropad, 0);
        for (size_t k = 0; k < frame_size; k++)
            fft_in[k] = window[k] * signal[k];
        // FFT
        fft_work = FFT::fftar_float(fftar_float_plan, block_size * zeropad, fft_in, fft_out, fft_work, fft_workN);
        std::copy(fft_out, fft_out + block_size * zeropad, out.begin());
        out[block_size * zeropad] = out[1];
        out[block_size * zeropad + 1] = 0;
        out[1] = 0;

        // subtract spectrum from audio spectrum
        for (size_t d = 0; d < block_size * zeropad; d += 2) {
            double re = out[d], im = out[d + 1];
            double sub_mag = sqrt(re * re + im * im);

            double mag = magnitude(audio_blocks[frame].noise.begin() + (long)d);
            if (mag > 0) {
                audio_blocks[frame].noise[d] /= (float)mag;
                audio_blocks[frame].noise[d + 1] /= (float)mag;
                mag -= sub_mag;
                if (mag < 0)
                    mag = 0;
                audio_blocks[frame].noise[d] *= (float)mag;
                audio_blocks[frame].noise[d + 1] *= (float)mag;
            }
        }

        if (killed("_subtract", frame & 7))
            return;
    }
    FFT::free_array_float(fft_in);
    FFT::free_array_float(fft_out);
    FFT::free_array_float(fft_work);
    FFT::cleanup(fftar_float_plan);
}

template <class AIter, class BIter> static double float_vector_delta(AIter ai, AIter aend, BIter bi) {
    double d = 0;
    while (ai != aend) {
        double dd = *ai++ - *bi++;
        d += dd * dd;
    }
    return d;
}

static void refine_sine_params_fast(EncoderBlock& audio_block, double mix_freq, const vector<float>& window) {
    const size_t frame_size = audio_block.debug_samples.size();

    AlignedArray<float, 16> sin_vec(frame_size);
    AlignedArray<float, 16> cos_vec(frame_size);
    AlignedArray<float, 16> sines(frame_size);
    AlignedArray<float, 16> all_sines(frame_size);

    vector<float> good_freqs;
    vector<float> good_mags;
    vector<float> good_phases;

    // figure out normalization for window
    double window_weight = 0;
    for (size_t i = 0; i < frame_size; i++)
        window_weight += window[i];

    for (size_t i = 0; i < audio_block.freqs.size(); i++) {
        VectorSinParams params;

        params.mix_freq = mix_freq;
        params.freq = audio_block.freqs[i];
        params.mag = audio_block.mags[i];
        params.phase = audio_block.phases[i];
        params.mode = VectorSinParams::ADD;

        fast_vector_sinf(params, &all_sines[0], &all_sines[frame_size]);
    }

    double max_mag;
    size_t partial = 0;
    do {
        max_mag = 0;
        // search biggest partial
        for (size_t i = 0; i < audio_block.freqs.size(); i++) {
            const double mag = audio_block.mags[i];

            if (mag > max_mag) {
                partial = i;
                max_mag = mag;
            }
        }
        // compute reconstruction of that partial
        if (max_mag > 0) {
            // remove partial, so we only do each partial once
            double f = audio_block.freqs[partial];

            audio_block.mags[partial] = 0;

            double phase;
            // determine "perfect" phase and magnitude instead of using interpolated fft phase
            double x_re = 0;
            double x_im = 0;

            VectorSinParams params;

            params.mix_freq = mix_freq;
            params.freq = f;
            params.mag = 1;
            params.phase = -((frame_size - 1) / 2.0) * f / mix_freq * 2.0 * M_PI;
            params.phase = normalize_phase(params.phase);
            params.mode = VectorSinParams::REPLACE;

            fast_vector_sincosf(params, &sin_vec[0], &sin_vec[frame_size], &cos_vec[0]);

            params.freq = f;
            params.mag = max_mag;
            params.phase = audio_block.phases[partial];
            params.mode = VectorSinParams::REPLACE;

            fast_vector_sinf(params, &sines[0], &sines[frame_size]);

            for (size_t n = 0; n < frame_size; n++) {
                double v = audio_block.debug_samples[n] - all_sines[n] + sines[n];
                v *= window[n];

                // multiply windowed signal with complex exp function from fourier transform:
                //
                //   v * exp (-j * x) = v * (cos (x) - j * sin (x))
                x_re += v * cos_vec[n];
                x_im -= v * sin_vec[n];
            }

            // correct influence of mirrored window (caused by negative frequency component)
            params.mix_freq = mix_freq;
            params.freq = 2 * f;
            params.mag = 1;
            params.phase = -((frame_size - 1) / 2.0) * (2 * f) / mix_freq * 2.0 * M_PI + 0.5 * M_PI;
            params.phase = normalize_phase(params.phase);
            params.mode = VectorSinParams::REPLACE;
            fast_vector_sinf(params, &cos_vec[0], &cos_vec[frame_size]);

            double w2omega = 0;
            for (size_t n = 0; n < frame_size; n++)
                w2omega += window[n] * cos_vec[n];

            x_re *= 2 / (window_weight + w2omega);
            x_im *= 2 / (window_weight - w2omega);

            // compute final magnitude & phase
            double magnitude = sqrt(x_re * x_re + x_im * x_im);
            phase = atan2(x_im, x_re) + 0.5 * M_PI;
            phase -= (frame_size - 1) / 2.0 / mix_freq * f * 2 * M_PI;
            phase = normalize_phase(phase);

            // restore partial => sines; keep params.freq & params.mix_freq
            params.freq = f;
            params.phase = phase;
            params.mag = magnitude;
            params.mode = VectorSinParams::ADD;
            fast_vector_sinf(params, &sines[0], &sines[frame_size]);

            // store refined freq, mag and phase
            good_freqs.push_back((float)f);
            good_mags.push_back((float)magnitude);
            good_phases.push_back((float)phase);
        }
    } while (max_mag > 0);

    audio_block.freqs = good_freqs;
    audio_block.mags = good_mags;
    audio_block.phases = good_phases;
}

static void remove_small_partials(EncoderBlock& audio_block) {
    /*
     * this function mainly serves to eliminate side peaks introduced by windowing
     * since these side peaks are typically much smaller than the main peak, we can
     * get rid of them by comparing peaks to the nearest peak, and removing them
     * if the nearest peak is much larger
     */
    vector<double> dbmags;
    for (vector<float>::iterator mi = audio_block.mags.begin(); mi != audio_block.mags.end(); mi++)
        dbmags.push_back(db_from_factor(*mi, -200));

    vector<bool> remove(dbmags.size());

    for (size_t i = 0; i < dbmags.size(); i++) {
        for (size_t j = 0; j < dbmags.size(); j++) {
            if (i != j) {
                double octaves = log(abs(audio_block.freqs[i] - audio_block.freqs[j])) / log(2);
                double mask = -30 - 15 * octaves; /* theoretical values -31 and -18 */
                if (dbmags[j] < dbmags[i] + mask)
                    ;
                // remove[j] = true;
            }
        }
    }

    vector<float> good_freqs;
    vector<float> good_mags;
    vector<float> good_phases;

    for (size_t i = 0; i < dbmags.size(); i++) {
        if (!remove[i]) {
            good_freqs.push_back(audio_block.freqs[i]);
            good_mags.push_back(audio_block.mags[i]);
            good_phases.push_back(audio_block.phases[i]);
        }
    }
    audio_block.freqs = good_freqs;
    audio_block.mags = good_mags;
    audio_block.phases = good_phases;
}

/**
 * This function reestimates the magnitudes and phases of the partials found
 * in the previous steps.
 */
void Encoder::optimize_partials(int optimization_level) {
    const double mix_freq = enc_params.mix_freq;

    for (uint64 frame = 0; frame < audio_blocks.size(); frame++) {
        if (optimization_level >= 1) // redo FFT estmates, only better
            refine_sine_params_fast(audio_blocks[frame], mix_freq, enc_params.window);

        remove_small_partials(audio_blocks[frame]);

        if (killed("_optimize"))
            return;
    }
}

static double mel_to_hz(double mel) {
    return 700 * (exp(mel / 1127.0) - 1);
}

static void approximate_noise_spectrum(double mix_freq, const vector<double>& spectrum, vector<double>& envelope,
                                       double norm) {
    size_t bands = envelope.size();
    size_t d = 0;
    for (size_t band = 0; band < envelope.size(); band++) {
        double mel_low = 30 + 4000.0 / bands * band;
        double mel_high = 30 + 4000.0 / bands * (band + 1);
        double hz_low = mel_to_hz(mel_low);
        double hz_high = mel_to_hz(mel_high);

        envelope[band] = 0;

        /* skip frequencies which are too low to be in lowest band */
        if (band == 0) {
            double f_hz = mix_freq / 2.0 * d / spectrum.size();
            while (f_hz < hz_low) {
                d += 2;
                f_hz = mix_freq / 2.0 * d / spectrum.size();
            }
        }
        double f_hz = mix_freq / 2.0 * d / spectrum.size();
        int n_values = 0;
        while (f_hz < hz_high) {
            if (d < spectrum.size()) {
                envelope[band] += (spectrum[d] * spectrum[d] + spectrum[d + 1] * spectrum[d + 1]);
                n_values++;
            }
            d += 2;
            f_hz = mix_freq / 2.0 * d / spectrum.size();
        }
        if (n_values > 0)
            envelope[band] = sqrt(envelope[band] / norm / n_values);
    }
}

/**
 * This function tries to approximate the residual by a spectral envelope
 * for a noise signal.
 */
void Encoder::approx_noise() {
    const size_t frame_size = enc_params.frame_size;
    const auto& window = enc_params.window;

    double sum_w2 = 0;
    for (size_t x = 0; x < frame_size; x++)
        sum_w2 += window[x] * window[x];

    // sum_w2 is the average influence of the window (w[x]^2), multiplied with frame_size
    const double norm = 0.5 * enc_params.mix_freq * sum_w2;

    for (uint64 frame = 0; frame < audio_blocks.size(); frame++) {
        vector<double> noise_envelope(32);
        vector<double> spectrum(audio_blocks[frame].noise.begin(), audio_blocks[frame].noise.end());

        /* A complex FFT would preserve the energy of the input signal exactly; the difference to
         * our (real) FFT is that every value in the complex spectrum occurs twice, once as "positive"
         * frequency, once as "negative" frequency - except for two spectrum values: the value
         * for frequency 0, and the value for frequency mix_freq / 2.
         *
         * To make this FFT energy preserving, we scale those values with a factor of sqrt (2) so
         * that their energy is twice as big (energy == squared value). Then we scale the whole
         * thing with a factor of 0.5, and we get an energy preserving transformation.
         */
        spectrum[0] /= sqrt(2);
        spectrum[spectrum.size() - 2] /= sqrt(2);

        approximate_noise_spectrum(enc_params.mix_freq, spectrum, noise_envelope, norm);
        audio_blocks[frame].noise.assign(noise_envelope.begin(), noise_envelope.end());

        if (killed("_noise", frame & 7))
            return;
    }
}

struct PartialData {
    float freq;
    float mag;
    float phase;
};

static bool pd_cmp(const PartialData& p1, const PartialData& p2) {
    return p1.freq < p2.freq;
}

void Encoder::sort_freqs() {
    for (uint64 frame = 0; frame < audio_blocks.size(); frame++) {
        // sort partials by frequency
        vector<PartialData> pvec;

        for (size_t p = 0; p < audio_blocks[frame].freqs.size(); p++) {
            PartialData pd;
            pd.freq = audio_blocks[frame].freqs[p];
            pd.mag = audio_blocks[frame].mags[p];
            pd.phase = audio_blocks[frame].phases[p];
            pvec.push_back(pd);
        }
        sort(pvec.begin(), pvec.end(), pd_cmp);

        // replace partial data with sorted partial data
        audio_blocks[frame].freqs.clear();
        audio_blocks[frame].mags.clear();
        audio_blocks[frame].phases.clear();

        for (vector<PartialData>::const_iterator pi = pvec.begin(); pi != pvec.end(); pi++) {
            // attack envelope computation produces some partials with mag = 0; we don't need to store these
            if (!approximatelyEqual(pi->mag, 0.0f)) {
                audio_blocks[frame].freqs.push_back(pi->freq);
                audio_blocks[frame].mags.push_back(pi->mag);
                audio_blocks[frame].phases.push_back(pi->phase);
            }
        }
    }
}

/**
 * This function calls all steps necessary for encoding in the right order.
 *
 * \param optimization_level determines if fast (0), medium (1), or very slow (2) algorithm is used
 */
bool Encoder::encode(const WavData& wav_data, int channel, int optimization_level, bool track_sines) {
    compute_stft(wav_data, channel);
    if (killed("stft"))
        return false;

    if (track_sines) {
        search_local_maxima();
        if (killed("maxima"))
            return false;

        link_partials();
        if (killed("link"))
            return false;

        validate_partials();
        if (killed("validate"))
            return false;

        optimize_partials(optimization_level);
        if (killed("optimize"))
            return false;

        spectral_subtract();
        if (killed("subtract"))
            return false;
    }
    approx_noise();
    if (killed("noise"))
        return false;

    sort_freqs();
    if (killed("sort"))
        return false;

    return true;
}

void Encoder::set_loop(Audio::LoopType loop_type_, int loop_start_, int loop_end_) {
    this->loop_type = loop_type_;
    this->loop_start = loop_start_;
    this->loop_end = loop_end_;
}

static void convert_freqs_mags_phases(const EncoderBlock& eblock, AudioBlock& ablock, const EncoderParams& enc_params) {
    const size_t frame_size = enc_params.frame_size;
    const double mix_freq = enc_params.mix_freq;
    const double fundamental_freq = enc_params.fundamental_freq;

    ablock.freqs.clear();
    ablock.mags.clear();
    ablock.phases.clear();

    for (size_t i = 0; i < eblock.freqs.size(); i++) {
        const uint16_t ifreq = sm_freq2ifreq(eblock.freqs[i] / fundamental_freq);
        const uint16_t imag = sm_factor2idb(eblock.mags[i]);

        double xphase = eblock.phases[i];

        // xphase is used (instead of phase) so that restoring the partial with the
        // quantized frequency will produce a minimal error at the center of the frame

        // compute unquantized phase at the center of the frame
        xphase += 2 * M_PI * eblock.freqs[i] / mix_freq * (frame_size - 1) / 2.0;

        // use quantized frequency to compute phase at the beginning of the frame
        xphase -= 2 * M_PI * sm_ifreq2freq(ifreq) * fundamental_freq / mix_freq * (frame_size - 1) / 2.0;

        // => normalize to interval [0..2*pi]
        xphase = normalize_phase(xphase);

        const uint16_t iphase = (uint16_t)sm_bound<int>(0, sm_round_positive(xphase / 2 / M_PI * 65536), 65535);

        // corner frequencies are most likely not part of the sound, but analysis
        // errors; so we don't save the smallest/largest possible freq
        if (ifreq != 0 && ifreq != 65535) {
            ablock.freqs.push_back(ifreq);
            ablock.mags.push_back(imag);

            if (enc_params.enable_phases)
                ablock.phases.push_back(iphase);
        }
    }
}

static void convert_noise(const vector<float>& noise, AudioBlock::Block& inoise) {
    inoise.resize(noise.size());

    for (size_t i = 0; i < noise.size(); i++)
        inoise[i] = sm_factor2idb(noise[i]);
}

/**
 * This function saves the data produced by the encoder to a SpectMorph file.
 */
Error Encoder::save(const string& filename) {
    std::unique_ptr<Audio> audio(save_as_audio());

    return audio->save(filename); // saving can fail
}

/**
 * This function saves the data produced by the encoder, returning a newly
 * allocated Audio object (caller must free this).
 */
Audio* Encoder::save_as_audio() {
    Audio* audio = new Audio();

    audio->fundamental_freq = (float)enc_params.fundamental_freq;
    audio->mix_freq = enc_params.mix_freq;
    audio->frame_size_ms = enc_params.frame_size_ms;
    audio->frame_step_ms = enc_params.frame_step_ms;
    audio->zero_values_at_start = (int)zero_values_at_start;
    audio->zeropad = (int)enc_params.zeropad;

    for (vector<EncoderBlock>::iterator ai = audio_blocks.begin(); ai != audio_blocks.end(); ai++) {
        AudioBlock block;
        convert_freqs_mags_phases(*ai, block, enc_params);
        convert_noise(ai->noise, block.noise);
        audio->contents.push_back(block);
    }
    audio->sample_count = (int)sample_count;
    audio->original_samples = original_samples;
    if (loop_start >= 0 && loop_end >= 0 && loop_type != Audio::LOOP_NONE) {
        audio->loop_type = loop_type;
        audio->loop_start = loop_start;
        audio->loop_end = loop_end;

        if (audio->loop_type == Audio::LOOP_TIME_FORWARD || audio->loop_type == Audio::LOOP_TIME_PING_PONG) {
            audio->loop_start += zero_values_at_start;
            audio->loop_end += zero_values_at_start;
        }
    }
    return audio;
}

/**
 * \mainpage SpectMorph Index Page
 *
 * \section intro_sec Introduction
 *
 * SpectMorph is a software which analyzes wav files and builds a frame based model of these wav files,
 * where each frame contains information about the spectrum. Each frame is represented as a sum of sine
 * waves, and a noise component. There are command line tools like smenc and smplay for encoding and
 * decoding SpectMorph models, which are documented in the manual pages. Technically, these tools are frontends
 * to the C++ classes in libspectmorph, which are documented here.
 *
 * \section enc_sec Encoding, loading and saving
 *
 * The encoder is implemented in SpectMorph::Encoder. It can be used to encode an audio file; the frames
 * (short snippets of the audio file, maybe 40 ms or so) are then available as vector containing
 * SpectMorph::AudioBlock objects.
 *
 * Many SpectMorph::AudioBlock objects are needed to represent a whole sound file, and the SpectMorph::Audio
 * class is used to store all parameters of an encoded file, along with the actual frames. Information like
 * the sampling rate (mix_freq) or the original note (fundamental freq) are stored in the SpectMorph::Audio
 * class. Functions for storing and loading SpectMorph::Audio objects exist in namespace SpectMorph::AudioFile.
 *
 * \section dec_sec Decoding
 *
 * The decoding process is frame oriented. For each SpectMorph::AudioBlock, a SpectMorph::Frame object needs
 * to be created. For each frame, a SpectMorph::SineDecoder and a SpectMorph::NoiseDecoder can be used to
 * reconstruct (something which sounds like) the original signal.
 */
