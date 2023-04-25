// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_AUDIO_HH
#define SPECTMORPH_AUDIO_HH

#include <vector>

#include "SMGenericIn.h"
#include "SMGenericOut.h"
#include "SMMath.h"
#include "SMUtils.h"

#define SPECTMORPH_BINARY_FILE_VERSION 140

namespace SpectMorph {

/**
 * \brief Block of audio data, encoded in SpectMorph parametric format
 *
 * This represents a single analysis frame, usually containing sine waves in freqs
 * and phases, and a noise envelope for everything that remained after subtracting
 * the sine waves.
 */
class AudioBlock {
  public:
    constexpr static int block_size = 1024;

    typedef std::vector<uint16_t> Block;
    Block noise;  //!< noise envelope, representing the original signal minus sine components
    Block freqs;  //!< frequencies of the sine components of this frame
    Block mags;   //!< magnitudes of the sine components
    Block phases; //!< phases of the sine components

    AudioBlock() {
        noise.reserve(block_size);
        freqs.reserve(block_size);
        mags.reserve(block_size);
        phases.reserve(block_size);
    }

    void sort_freqs();
    double estimate_fundamental(int n_partials = 1, double* mag = nullptr) const;

    double freqs_f(size_t i) const {
        return sm_ifreq2freq(freqs[i]);
    }

    double mags_f(size_t i) const {
        return sm_idb2factor(mags[i]);
    }

    double phases_f(size_t i) const {
        const double factor = 2.0 * M_PI / 65536.0;
        return phases[i] * factor;
    }

    double noise_f(size_t i) const {
        return sm_idb2factor(noise[i]);
    }
};

/**
 * \brief Audio sample containing many blocks
 *
 * This class contains the information the SpectMorph::Encoder creates for a wav file. The
 * time dependant parameters are stored in contents, as a vector of audio frames; the
 * parameters that are the same for all frames are stored in this class.
 */
class Audio {
    SPECTMORPH_CLASS_NON_COPYABLE(Audio);

  public:
    Audio();
    ~Audio();

    enum LoopType {
        LOOP_NONE = 0,
        LOOP_FRAME_FORWARD,
        LOOP_FRAME_PING_PONG,
        LOOP_TIME_FORWARD,
        LOOP_TIME_PING_PONG,
    };

    float fundamental_freq = 0;          //!< fundamental frequency (note which was encoded), or 0 if not available
    float mix_freq = 0;                  //!< mix freq (sampling rate) of the original audio data
    float frame_size_ms = 0;             //!< length of each audio frame in milliseconds
    float frame_step_ms = 0;             //!< stepping of the audio frames in milliseconds
    int zeropad = 0;                     //!< FFT zeropadding used during analysis
    LoopType loop_type = LOOP_NONE;      //!< type of loop to be used during sustain phase of playback
    int start = 0;                       //!< playback start position
    int end = 0;                         //!< playback end position
    int loop_start = 0;                  //!< loop point to be used during sustain phase of playback
    int loop_end = 0;                    //!< loop point to be used during sustain phase of playback
    int zero_values_at_start = 0;        //!< number of zero values added by encoder (strip during decoding)
    int sample_count = 0;                //!< number of samples encoded (including zero_values_at_start)
    std::vector<float> original_samples; //!< original time domain signal as samples (debugging only)
    float original_samples_norm_db = 0;  //!< normalization factor to be applied to original samples
    int original_start = 0;              //!< original sample start in samples
    int original_end = 0;                //!< original sample end in samples
    int original_loop_start = 0;         //!< original loop start in samples
    int original_loop_end = 0;           //!< original loop end in samples
    std::vector<AudioBlock> contents;    //!< the actual frame data

    Error load(const std::string& filename);
    Error load(SpectMorph::GenericIn* file);
    Error save(const std::string& filename) const;
    Error save(SpectMorph::GenericOut* file) const;

    static bool loop_type_to_string(LoopType loop_type, std::string& s);
    static bool string_to_loop_type(const std::string& s, LoopType& loop_type);
};

} // namespace SpectMorph

#endif
