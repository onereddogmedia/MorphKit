// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_FFT_HH
#define SPECTMORPH_FFT_HH

#include "pffft.h"
#include <map>
#include <sys/types.h>

namespace SpectMorph {

namespace FFT {

float* new_array_float(size_t N);
void free_array_float(float* f);

typedef std::map<size_t, PFFFT_Setup*> PlanMap;

float* fftar_float(PlanMap& plan_map, size_t N, float* in, float* out, float* work, size_t& workN);
void fftsr_destructive_float(PlanMap& plan_map, size_t N, float* in, float* out, float* work);
void cleanup(PlanMap& plan_map);

} // namespace FFT

} // namespace SpectMorph

#endif
