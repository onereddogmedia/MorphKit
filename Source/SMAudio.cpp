// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMAudio.h"
#include "SMInFile.h"
#include "SMLeakDebugger.h"
#include "SMMMapIn.h"
#include "SMMemOut.h"
#include "SMOutFile.h"
#include "SMStdioOut.h"
#include <assert.h>
#include <errno.h>
#include <fcntl.h>
#include <stdio.h>

using std::string;
using std::vector;

using namespace SpectMorph;

static LeakDebugger leak_debugger("SpectMorph::Audio");

/**
 * This function loads a SM-File.
 *
 * \param filename the name of the SM-File to be loaded
 * \returns a SpectMorph::Error indicating whether loading was successful
 */
Error SpectMorph::Audio::load(const string& filename) {
    GenericIn* file = GenericIn::open(filename);
    if (!file)
        return Error::Code::FILE_NOT_FOUND;

    Error result = load(file);
    delete file;

    return result;
}

Error SpectMorph::Audio::load(GenericIn* file) {
    SpectMorph::AudioBlock* audio_block = nullptr;

    InFile ifile(file);

    string section;
    size_t contents_pos = 0;

    if (!ifile.open_ok())
        return Error::Code::FILE_NOT_FOUND;

    if (ifile.file_type() != "SpectMorph::Audio")
        return Error::Code::FORMAT_INVALID;

    if (ifile.file_version() != SPECTMORPH_BINARY_FILE_VERSION)
        return Error::Code::FORMAT_INVALID;

    while (ifile.event() != InFile::END_OF_FILE) {
        if (ifile.event() == InFile::BEGIN_SECTION) {
            assert(section == "");
            section = ifile.event_name();

            if (section == "frame") {
                assert(audio_block == nullptr);
                assert(contents_pos < contents.size());

                audio_block = &contents[contents_pos];
            }
        } else if (ifile.event() == InFile::END_SECTION) {
            if (section == "frame") {
                assert(audio_block);

                contents_pos++;
                audio_block = nullptr;
            }

            assert(section != "");
            section = "";
        } else if (ifile.event() == InFile::INT) {
            if (section == "header") {
                if (ifile.event_name() == "zeropad")
                    zeropad = ifile.event_int();
                else if (ifile.event_name() == "start")
                    start = ifile.event_int();
                else if (ifile.event_name() == "end")
                    end = ifile.event_int();
                else if (ifile.event_name() == "loop_start")
                    loop_start = ifile.event_int();
                else if (ifile.event_name() == "loop_end")
                    loop_end = ifile.event_int();
                else if (ifile.event_name() == "loop_type")
                    loop_type = static_cast<LoopType>(ifile.event_int());
                else if (ifile.event_name() == "zero_values_at_start")
                    zero_values_at_start = ifile.event_int();
                else if (ifile.event_name() == "sample_count")
                    sample_count = ifile.event_int();
                else if (ifile.event_name() == "original_start")
                    original_start = ifile.event_int();
                else if (ifile.event_name() == "original_end")
                    original_end = ifile.event_int();
                else if (ifile.event_name() == "original_loop_start")
                    original_loop_start = ifile.event_int();
                else if (ifile.event_name() == "original_loop_end")
                    original_loop_end = ifile.event_int();
                else if (ifile.event_name() == "frame_count") {
                    int frame_count = ifile.event_int();
                    contents.clear();
                    contents.resize((size_t)frame_count);
                    contents_pos = 0;
                } else {
                    printf("unhandled int %s %s\n", section.c_str(), ifile.event_name().c_str());
                }
            } else {
                assert(false);
            }
        } else if (ifile.event() == InFile::FLOAT) {
            if (section == "header") {
                if (ifile.event_name() == "mix_freq")
                    mix_freq = ifile.event_float();
                else if (ifile.event_name() == "frame_size_ms")
                    frame_size_ms = ifile.event_float();
                else if (ifile.event_name() == "frame_step_ms")
                    frame_step_ms = ifile.event_float();
                else if (ifile.event_name() == "fundamental_freq")
                    fundamental_freq = ifile.event_float();
                else if (ifile.event_name() == "original_samples_norm_db")
                    original_samples_norm_db = ifile.event_float();
                else
                    printf("unhandled float %s  %s\n", section.c_str(), ifile.event_name().c_str());
            } else {
                assert(false);
            }
        } else if (ifile.event() == InFile::FLOAT_BLOCK) {
            const vector<float>& fb = ifile.event_float_block();

            if (section == "header") {
                if (ifile.event_name() == "original_samples") {
                    original_samples = fb;
                } else {
                    printf("unhandled float block %s  %s\n", section.c_str(), ifile.event_name().c_str());
                }
            }
        } else if (ifile.event() == InFile::UINT16_BLOCK) {
            const vector<uint16_t>& ib = ifile.event_uint16_block();
            if (ifile.event_name() == "freqs") {
                std::copy(ib.begin(), ib.end(), std::back_inserter(audio_block->freqs));

                // ensure that freqs are sorted (we need that for LiveDecoder)
                int old_freq = -1;

                for (size_t i = 0; i < ib.size(); i++) {
                    if (ib[i] < old_freq) {
                        printf("frequency data is not sorted, can't play file\n");
                        return Error::Code::PARSE_ERROR;
                    }
                    old_freq = ib[i];
                }
            } else if (ifile.event_name() == "mags") {
                std::copy(ib.begin(), ib.end(), std::back_inserter(audio_block->mags));
            } else if (ifile.event_name() == "phases") {
                std::copy(ib.begin(), ib.end(), std::back_inserter(audio_block->phases));
            } else if (ifile.event_name() == "noise") {
                std::copy(ib.begin(), ib.end(), std::back_inserter(audio_block->noise));
            } else {
                printf("unhandled int16 block %s %s\n", section.c_str(), ifile.event_name().c_str());
                assert(false);
            }
        } else if (ifile.event() == InFile::READ_ERROR) {
            return Error::Code::PARSE_ERROR;
        } else {
            return Error::Code::PARSE_ERROR;
        }
        ifile.next_event();
    }
    return Error::Code::NONE;
}

SpectMorph::Audio::Audio() {
    leak_debugger.add(this);
}

Audio::~Audio() {
    leak_debugger.del(this);
}

/**
 * This function saves a SM-File.
 *
 * \param filename the name of the SM-File to be written
 * \returns a SpectMorph::Error indicating saving loading was successful
 */
Error SpectMorph::Audio::save(const string& filename) const {
    GenericOut* out = StdioOut::open(filename);
    if (!out) {
        fprintf(stderr, "error: can't open output file '%s'.\n", filename.c_str());
        return Error::Code::FILE_NOT_FOUND;
    }
    Error result = save(out);
    delete out; // close file

    return result;
}

Error SpectMorph::Audio::save(GenericOut* file) const {
    OutFile of(file, "SpectMorph::Audio", SPECTMORPH_BINARY_FILE_VERSION);
    assert(of.open_ok());

    of.begin_section("header");
    of.write_float("mix_freq", mix_freq);
    of.write_float("frame_size_ms", frame_size_ms);
    of.write_float("frame_step_ms", frame_step_ms);
    of.write_float("fundamental_freq", fundamental_freq);
    of.write_float("original_samples_norm_db", original_samples_norm_db);
    of.write_int("zeropad", zeropad);
    of.write_int("start", start);
    of.write_int("end", end);
    of.write_int("loop_type", loop_type);
    of.write_int("loop_start", loop_start);
    of.write_int("loop_end", loop_end);
    of.write_int("zero_values_at_start", zero_values_at_start);
    of.write_int("frame_count", (int)contents.size());
    of.write_int("sample_count", sample_count);
    of.write_int("original_start", original_start);
    of.write_int("original_end", original_end);
    of.write_int("original_loop_start", original_loop_start);
    of.write_int("original_loop_end", original_loop_end);
    of.write_float_block("original_samples", original_samples);
    of.end_section();

    for (size_t i = 0; i < contents.size(); i++) {
        // ensure that freqs are sorted (we need that for LiveDecoder)
        int old_freq = -1;
        for (size_t f = 0; f < contents[i].freqs.size(); f++) {
            assert(contents[i].freqs[f] >= old_freq);
            old_freq = contents[i].freqs[f];
        }

        of.begin_section("frame");
        of.write_uint16_block("noise", contents[i].noise);
        of.write_uint16_block("freqs", contents[i].freqs);
        of.write_uint16_block("mags", contents[i].mags);
        of.write_uint16_block("phases", contents[i].phases);
        of.end_section();
    }
    return Error::Code::NONE;
}

bool Audio::loop_type_to_string(LoopType loop_type, string& s) {
    switch (loop_type) {
        case LOOP_NONE: {
            s = "loop-none";
            break;
        }
        case LOOP_FRAME_FORWARD: {
            s = "loop-frame-forward";
            break;
        }
        case LOOP_FRAME_PING_PONG: {
            s = "loop-frame-ping-pong";
            break;
        }
        case LOOP_TIME_FORWARD: {
            s = "loop-time-forward";
            break;
        }
        case LOOP_TIME_PING_PONG: {
            s = "loop-time-ping-pong";
            break;
        }
        default: {
            return false; // unknown loop type
        }
    }
    return true;
}

bool Audio::string_to_loop_type(const string& s, LoopType& loop_type) {
    if (s == "loop-none") {
        loop_type = LOOP_NONE;
    } else if (s == "loop-frame-forward") {
        loop_type = LOOP_FRAME_FORWARD;
    } else if (s == "loop-frame-ping-pong") {
        loop_type = LOOP_FRAME_PING_PONG;
    } else if (s == "loop-time-forward") {
        loop_type = LOOP_TIME_FORWARD;
    } else if (s == "loop-time-ping-pong") {
        loop_type = LOOP_TIME_PING_PONG;
    } else {
        return false; // unknown loop type
    }
    return true;
}

namespace {

struct PartialData {
    uint16_t freq;
    uint16_t mag;
};

static bool pd_cmp(const PartialData& p1, const PartialData& p2) {
    return p1.freq < p2.freq;
}

} // namespace

void AudioBlock::sort_freqs() {
    // sorting is required for morphing generated blocks only, which have no phase information
    g_return_if_fail(phases.empty());

    // sort partials by frequency
    const size_t N = freqs.size();
    PartialData pvec[N + AVOID_ARRAY_UB];

    for (size_t p = 0; p < N; p++) {
        pvec[p].freq = freqs[p];
        pvec[p].mag = mags[p];
    }
    std::sort(pvec, pvec + N, pd_cmp);

    // replace partial data with sorted partial data
    for (size_t p = 0; p < N; p++) {
        freqs[p] = pvec[p].freq;
        mags[p] = pvec[p].mag;
    }
}

double AudioBlock::estimate_fundamental(int n_partials, double* mag) const {
    g_return_val_if_fail(n_partials >= 1 && n_partials <= 3, 1.0);

    double est_freq = 0, est_mag = 0;

    auto update_estimate = [&](int n, double freq_min, double freq_max) {
        if (n > n_partials)
            return;

        double best_freq = 0, best_mag = 0;

        for (size_t p = 0; p < mags.size(); p++) {
            if (freqs_f(p) > freq_min && freqs_f(p) < freq_max && mags_f(p) > best_mag) {
                best_mag = mags_f(p);
                best_freq = freqs_f(p) / n;
            }
        }
        if (best_mag > 0) {
            est_mag += best_mag;
            est_freq += best_freq * best_mag;
        }
    };

    update_estimate(1, 0.8, 1.25);
    update_estimate(2, 1.5, 2.5);
    update_estimate(3, 2.5, 3.5);

    if (mag)
        *mag = est_mag;

    if (est_mag > 0)
        return est_freq / est_mag;
    else
        return 1;
}
