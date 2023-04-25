// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_BLOCK_UTILS_HH
#define SPECTMORPH_BLOCK_UTILS_HH

#include "SMUtils.h"

namespace SpectMorph {

/* Block utils */

class Block {
  public:
    static void mul(uint n_values, float* ovalues, const float* ivalues);
    static void add(uint n_values, float* ovalues, const float* ivalues);
    static void zero(uint n_values, float* ovalues);
};

} // namespace SpectMorph

#endif /* SPECTMORPH_BLOCK_UTILS_HH */
