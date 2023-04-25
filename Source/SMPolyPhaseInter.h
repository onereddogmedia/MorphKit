// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_POLY_PHASE_INTER_HH
#define SPECTMORPH_POLY_PHASE_INTER_HH

#include "fixed_capacity_vector.h"
#include <sys/types.h>
#include <vector>

namespace SpectMorph {

class PolyPhaseInter {
  public:
    PolyPhaseInter();
    ~PolyPhaseInter() = default;

    typedef std::experimental::fixed_capacity_vector<float, 32768> Buffer;

    double get_sample(const Buffer& signal, double pos);
    double get_sample_no_check(const Buffer& signal, double pos);

    size_t get_min_padding();

  private:
    std::vector<float> x;
};

} // namespace SpectMorph

#endif
