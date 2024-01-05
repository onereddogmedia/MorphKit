// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_POLY_PHASE_INTER_HH
#define SPECTMORPH_POLY_PHASE_INTER_HH

#include <sys/types.h>
#include <vector>

namespace SpectMorph {

class PolyPhaseInter {
    PolyPhaseInter();
    ~PolyPhaseInter() {
    }

    std::vector<float> x;

  public:
    static PolyPhaseInter* the();

    double get_sample(const std::vector<float>& signal, double pos);
    double get_sample_no_check(const float* signal, double pos);

    size_t get_min_padding();
};

} // namespace SpectMorph

#endif
