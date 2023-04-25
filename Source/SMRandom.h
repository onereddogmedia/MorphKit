// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_RANDOM_HH
#define SPECTMORPH_RANDOM_HH

#include <stdio.h>
#include <stdlib.h>

#include "SMPcg32rng.h"

namespace SpectMorph {

class Random {
    Pcg32Rng rand_gen;

  public:
    Random();

    void set_seed(uint32_t seed);

    inline double random_double_range(double begin, double end) {
        const uint32_t rand_max = 0xffffffff; // Pcg32Rng output: complete 32-bit values
        const uint32_t r = random_uint32();
        const double scale = 1.0 / (double(rand_max) + 1.0);

        return r * scale * (end - begin) + begin;
    }
    inline uint32_t random_uint32() {
        return rand_gen.random();
    }
    inline void random_block(size_t n_values, uint32_t* values) {
        while (n_values--)
            *values++ = random_uint32();
    }
};

class RandomRange {
  public:
    static void Update() {
        // Galois LFSR with feedback polynomial = x^16 + x^14 + x^13 + x^11.
        // Period: 65535.
        rng_state_ = (rng_state_ >> 1) ^ (-(rng_state_ & 1) & 0xb400);
    }

    static inline uint16_t state() {
        return rng_state_;
    }

    static inline void Seed(uint16_t seed) {
        rng_state_ = seed;
    }

    static inline uint8_t state_msb() {
        return static_cast<uint8_t>(rng_state_ >> 8);
    }

    static inline uint8_t GetByte() {
        Update();
        return state_msb();
    }

    static inline uint16_t GetWord() {
        Update();
        return state();
    }

    static inline float GetFloat() {
        return static_cast<float>(GetWord()) / 4294967296.0f;
    }

    static inline float GetRange(float begin, float end) {
        float r = GetFloat();
        return r * end - (r - 1) * begin;
    }

  private:
    static uint16_t rng_state_;
};

} // namespace SpectMorph
#endif
