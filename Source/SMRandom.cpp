// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMRandom.h"

#include <string.h>

#include "glib.h"

using SpectMorph::Random;

uint16_t SpectMorph::RandomRange::rng_state_ = 0x21;

Random::Random() {
    set_seed(3320);
}

void Random::set_seed(uint32_t seed) {
    /* multiply seed with different primes in order to set the generator state to
     * a pseudo-random position derived from seed
     */
    const uint64_t prime1 = 3126986573;
    const uint64_t prime2 = 4151919467;

    rand_gen.seed(prime1 * seed, prime2 * seed);
}
