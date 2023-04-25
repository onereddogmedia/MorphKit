// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMLiveDecoder.h"
#include "SMLeakDebugger.h"
#include "SMMath.h"
#include "SMUtils.h"

#include <assert.h>
#include <stdio.h>

using namespace SpectMorph;

using std::max;
using std::min;
using std::vector;

static LeakDebugger leak_debugger("SpectMorph::LiveDecoder");

#define ANTIALIAS_FILTER_TABLE_SIZE 256

static vector<float> antialias_filter_table;

static void init_aa_filter() {
    if (antialias_filter_table.empty()) {
        antialias_filter_table.resize(ANTIALIAS_FILTER_TABLE_SIZE);

        const double db_at_nyquist = -60;

        for (size_t i = 0; i < antialias_filter_table.size(); i++)
            antialias_filter_table[i] = (float)db_to_factor(double(i) / ANTIALIAS_FILTER_TABLE_SIZE * db_at_nyquist);
    }
}

static inline double fmatch(double f1, double f2) {
    return f2 < (f1 * 1.05) && f2 > (f1 * 0.95);
}

LiveDecoder::LiveDecoder() : audio(nullptr), ifft_synth(nullptr), source(nullptr), sse_samples(nullptr) {
    init_aa_filter();
    pp_inter = new PolyPhaseInter();
    leak_debugger.add(this);
}

LiveDecoder::LiveDecoder(LiveDecoderSource* source_) : LiveDecoder() {
    this->source = source_;
}

LiveDecoder::~LiveDecoder() {
    if (ifft_synth) {
        delete ifft_synth;
        ifft_synth = nullptr;
    }
    if (sse_samples) {
        delete sse_samples;
        sse_samples = nullptr;
    }
    if (pp_inter) {
        delete pp_inter;
        pp_inter = nullptr;
    }
    leak_debugger.del(this);
}

static size_t preferred_block_size(double mix_freq) {
    size_t bs = 1;

    while (bs * 2 / mix_freq < 0.040) /* block size should not exceed 40ms (and be power of 2) */
        bs *= 2;

    return bs;
}

void LiveDecoder::prepareToPlay(float mix_freq) {
    current_mix_freq = mix_freq;

    if (source) {
        source->prepareToPlay(mix_freq);
    }
    block_size = preferred_block_size(mix_freq);

    if (ifft_synth)
        delete ifft_synth;
    ifft_synth = new IFFTSynth(block_size, mix_freq, IFFTSynth::WIN_HANNING);

    if (sse_samples)
        delete sse_samples;
    sse_samples = new AlignedArray<float, 16>(block_size);
}

void LiveDecoder::retrigger(int channel, float freq, int midi_velocity) {
    Audio* best_audio = nullptr;

    if (source) {
        source->retrigger(channel, freq, midi_velocity);
        best_audio = source->audio();
    }
    audio = best_audio;

    if (best_audio) {
        frame_size = (size_t)(audio->frame_size_ms * current_mix_freq / 1000);
        frame_step = (size_t)(audio->frame_step_ms * current_mix_freq / 1000);
        zero_values_at_start_scaled = (size_t)(audio->zero_values_at_start * current_mix_freq / audio->mix_freq);
        loop_start_scaled = (size_t)(audio->loop_start * current_mix_freq / audio->mix_freq);
        loop_end_scaled = (size_t)(audio->loop_end * current_mix_freq / audio->mix_freq);
        loop_point = (get_loop_type() == Audio::LOOP_NONE) ? -1 : audio->loop_start;

        have_samples = 0;
        pos = 0;
        frame_idx = 0;
        env_pos = 0;
        original_samples_norm_factor = db_to_factor(audio->original_samples_norm_db);

        // reset partial state vectors
        pstate[0].clear();
        pstate[1].clear();
        last_pstate = &pstate[0];

        // setup portamento state
        assert(PortamentoState::DELTA >= pp_inter->get_min_padding());

        portamento_state.pos = PortamentoState::DELTA;
        portamento_state.buffer.resize(PortamentoState::DELTA);
    }
    current_freq = freq;
}

size_t LiveDecoder::compute_loop_frame_index(size_t frame_idx, Audio* audio) {
    if (int(frame_idx) > audio->loop_start) {
        g_return_val_if_fail(audio->loop_end >= audio->loop_start, frame_idx);

        if (audio->loop_type == Audio::LOOP_FRAME_FORWARD) {
            size_t loop_len = (size_t)(audio->loop_end + 1 - audio->loop_start);
            frame_idx = (size_t)audio->loop_start + (frame_idx - (size_t)audio->loop_start) % loop_len;
        } else if (audio->loop_type == Audio::LOOP_FRAME_PING_PONG) {
            size_t loop_len = (size_t)(audio->loop_end - audio->loop_start);
            if (loop_len > 0) {
                size_t ping_pong_len = loop_len * 2;
                size_t ping_pong_pos = (frame_idx - (size_t)audio->loop_start) % ping_pong_len;

                if (ping_pong_pos < loop_len) // ping part of the ping-pong loop (forward)
                {
                    frame_idx = (size_t)audio->loop_start + ping_pong_pos;
                } else // pong part of the ping-pong loop (backward)
                {
                    frame_idx = (size_t)audio->loop_end - (ping_pong_pos - loop_len);
                }
            } else {
                frame_idx = (size_t)audio->loop_start;
            }
        }
    }
    return frame_idx;
}

void LiveDecoder::process_internal(size_t n_values, float* audio_out, float portamento_stretch) {
    assert(audio); // need selected (triggered) audio to use this function

    const double portamento_env_step = 1 / portamento_stretch;
    unsigned int i = 0;
    while (i < n_values) {
        if (have_samples == 0) {
            double want_freq = current_freq;

            std::copy(&(*sse_samples)[block_size / 2], &(*sse_samples)[block_size], &(*sse_samples)[0]);
            zero_float_block(block_size / 2, &(*sse_samples)[block_size / 2]);

            if (get_loop_type() == Audio::LOOP_TIME_FORWARD) {
                size_t xenv_pos = (size_t)env_pos;

                if (xenv_pos > loop_start_scaled) {
                    xenv_pos = (xenv_pos - loop_start_scaled) % (loop_end_scaled - loop_start_scaled);
                    xenv_pos += loop_start_scaled;
                }
                frame_idx = xenv_pos / frame_step;
            } else if (get_loop_type() == Audio::LOOP_FRAME_FORWARD || get_loop_type() == Audio::LOOP_FRAME_PING_PONG) {
                frame_idx = compute_loop_frame_index((size_t)(env_pos / frame_step), audio);
            } else {
                frame_idx = (size_t)(env_pos / frame_step);
                if (loop_point != -1 && frame_idx > size_t(loop_point)) /* if in loop mode: loop current frame */
                    frame_idx = (size_t)loop_point;
            }

            AudioBlock* audio_block_ptr = nullptr;
            if (source) {
                audio_block_ptr = source->audio_block(frame_idx);
            } else if (frame_idx < audio->contents.size()) {
                audio_block_ptr = &audio->contents[frame_idx];
            }
            if (audio_block_ptr) {
                const AudioBlock& audio_block = *audio_block_ptr;

                assert(audio_block.freqs.size() == audio_block.mags.size());

                ifft_synth->clear_partials();

                // point n_pstate to pstate[0] and pstate[1] alternately (one holds points to last state and the other
                // points to new state)
                bool lps_zero = (last_pstate == &pstate[0]);
                auto& new_pstate = lps_zero ? pstate[1] : pstate[0];
                const auto& old_pstate = lps_zero ? pstate[0] : pstate[1];

                new_pstate.clear(); // clear old partial state

                const double phase_factor = block_size * M_PI / current_mix_freq;
                const double filter_fact =
                    18000.0 / current_mix_freq; // filter at 18 kHz (higher mix freq => higher filter)
                const double filter_min_freq = filter_fact * current_mix_freq;

                size_t old_partial = 0;
                for (size_t partial = 0; partial < audio_block.freqs.size(); partial++) {
                    const double freq = audio_block.freqs_f(partial) * want_freq;

                    // anti alias filter:
                    double mag = audio_block.mags_f(partial);
                    double phase = 0; // atan2 (smag, cmag); FIXME: Does initial phase matter? I think not.

                    // portamento:
                    //  - portamento_stretch > 1 means we read out faster
                    //  => this means the aliasing starts at lower frequencies
                    const double portamento_freq = freq * max(portamento_stretch, 1.0f);
                    if (portamento_freq > filter_min_freq) {
                        double norm_freq = portamento_freq / current_mix_freq;
                        if (norm_freq > 0.5) {
                            // above nyquist freq -> since partials are sorted, there is nothing more to do for this
                            // frame
                            break;
                        } else {
                            // between filter_fact and 0.5 (db linear filter)
                            int index = sm_round_positive(ANTIALIAS_FILTER_TABLE_SIZE * (norm_freq - filter_fact) /
                                                          (0.5 - filter_fact));
                            if (index >= 0) {
                                if (index < ANTIALIAS_FILTER_TABLE_SIZE)
                                    mag *= (double)antialias_filter_table[(size_t)index];
                                else
                                    mag = 0;
                            }
                        }
                    }

                    /*
                     * increment old_partial as long as there is a better candidate (closer to freq)
                     */
                    bool freq_match = false;
                    if (!old_pstate.empty()) {
                        double best_fdiff = fabs(old_pstate[old_partial].freq - freq);

                        while ((old_partial + 1) < old_pstate.size()) {
                            double fdiff = fabs(old_pstate[old_partial + 1].freq - freq);
                            if (fdiff < best_fdiff) {
                                old_partial++;
                                best_fdiff = fdiff;
                            } else {
                                break;
                            }
                        }
                        const double lfreq = old_pstate[old_partial].freq;
                        freq_match = fmatch(lfreq, freq) > 0;
                    }

                    if (freq_match) {
                        // matching freq -> compute new phase
                        const double lfreq = old_pstate[old_partial].freq;
                        const double lphase = old_pstate[old_partial].phase;

                        phase = fmod(lphase + lfreq * phase_factor, 2 * M_PI);
                    }
                    ifft_synth->render_partial(freq, mag, phase);

                    PartialState ps;
                    ps.freq = (float)freq;
                    ps.phase = (float)phase;
                    new_pstate.push_back(ps);
                }
                last_pstate = &new_pstate;

                float* samples = &(*sse_samples)[0];
                ifft_synth->get_samples(samples, IFFTSynth::ADD);
            }
            pos = 0;
            have_samples = block_size / 2;
        }

        if (env_pos >= zero_values_at_start_scaled) {
            size_t can_copy = min(have_samples, n_values - i);

            memcpy(audio_out + i, &(*sse_samples)[pos], sizeof(float) * can_copy);
            i += can_copy;
            pos += can_copy;
            env_pos += can_copy * portamento_env_step;
            have_samples -= can_copy;
        } else {
            // skip sample
            pos++;
            env_pos += portamento_env_step;
            have_samples--;
        }
    }
}

static bool portamento_check(size_t n_values, const float* freq_in, float current_freq) {
    /* if any value in freq_in differs from current_freq => need portamento */
    for (size_t i = 0; i < n_values; i++) {
        if (fabs(freq_in[i] / current_freq - 1) > 0.0001) // very small frequency difference (less than one cent)
            return true;
    }

    /* all freq_in[i] are (approximately) current_freq => no portamento */
    return false;
}

void LiveDecoder::portamento_grow(double end_pos, float portamento_stretch) {
    /* produce input samples until current_pos */
    const int TODO = int(end_pos) + PortamentoState::DELTA - int(portamento_state.buffer.size());
    if (TODO > 0) {
        const size_t START = portamento_state.buffer.size();

        portamento_state.buffer.resize(portamento_state.buffer.size() + (size_t)TODO);
        process_internal((size_t)TODO, &portamento_state.buffer[START], portamento_stretch);
    }
    portamento_state.pos = end_pos;
}

void LiveDecoder::portamento_shrink() {
    auto& buffer = portamento_state.buffer;

    /* avoid infinite state */
    if (buffer.size() > 256) {
        const int shrink_buffer = (int)buffer.size() - 2 * PortamentoState::DELTA; // only keep 2 * DELTA samples

        buffer.erase(buffer.begin(), buffer.begin() + shrink_buffer);
        portamento_state.pos -= shrink_buffer;
    }
}

void LiveDecoder::process_portamento(size_t n_values, const float* freq_in, float* audio_out) {
    assert(audio); // need selected (triggered) audio to use this function

    const double start_pos = portamento_state.pos;
    const auto& buffer = portamento_state.buffer;

    if (!portamento_state.active) {
        if (freq_in && portamento_check(n_values, freq_in, current_freq))
            portamento_state.active = true;
    }
    if (portamento_state.active) {
        const size_t max_n_values = (const size_t)(kMaxSampleRate * 0.010f);

        float fake_freq_in[max_n_values];
        if (!freq_in) {
            std::fill(fake_freq_in, fake_freq_in + n_values, current_freq);
            freq_in = fake_freq_in;
        }

        double pos_[max_n_values], end_pos = start_pos, current_step = 1;

        for (size_t i = 0; i < n_values; i++) {
            pos_[i] = end_pos;

            current_step = freq_in[i] / current_freq;
            end_pos += current_step;
        }
        portamento_grow(end_pos, (float)current_step);

        /* interpolate from buffer (portamento) */
        for (size_t i = 0; i < n_values; i++)
            audio_out[i] = (float)pp_inter->get_sample_no_check(buffer, pos_[i]);
    } else {
        /* no portamento: just compute & copy values */
        portamento_grow(start_pos + n_values, 1);

        const float* start = &buffer[(size_t)sm_round_positive(start_pos)];
        std::copy(start, start + n_values, audio_out);
    }
    portamento_shrink();
}

void LiveDecoder::process(size_t n_values, const float* freq_in, float* audio_out) {
    if (!audio) // nothing loaded
    {
        std::fill(audio_out, audio_out + n_values, 0);
        return;
    }
    /* ensure that time_offset_ms() is only called during live decoder process */
    start_env_pos = env_pos;
    /*
     * split processing into small blocks, max 10 ms
     *  -> limit n_values to keep portamento stretch settings up-to-date
     */
    const size_t max_n_values = (const size_t)(current_mix_freq * 0.010f);

    while (n_values > 0 && sse_samples && ifft_synth) {
        size_t todo_values = min(n_values, max_n_values);

        process_portamento(todo_values, freq_in, audio_out);

        if (freq_in)
            freq_in += todo_values;

        audio_out += todo_values;
        n_values -= todo_values;
    }
}

Audio::LoopType LiveDecoder::get_loop_type() {
    assert(audio);

    return audio->loop_type;
}

double LiveDecoder::current_pos() const {
    if (!audio)
        return -1;

    return frame_idx * audio->frame_step_ms - audio->zero_values_at_start * 1000.0 / audio->mix_freq;
}

double LiveDecoder::fundamental_note() const {
    if (!audio)
        return -1;

    return sm_freq_to_note(audio->fundamental_freq);
}
