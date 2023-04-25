// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMFft.h"
#include "SMUtils.h"
#include "SMRandom.h"
#include <algorithm>

#include "glib.h"

using namespace SpectMorph;

static PFFFT_Setup* read_plan_map_threadsafe(FFT::PlanMap& plan_map, size_t N) {
    /* std::map access is not threadsafe */
    if (plan_map[N] == nullptr) {
        plan_map[N] = pffft_new_setup((int)N, PFFFT_REAL);
    }
    return plan_map[N];
}

float* FFT::new_array_float(size_t N) {
    const size_t N_2 = N + 2; /* extra space for r2c extra complex output */

    float* result = (float*)pffft_aligned_malloc(sizeof(float) * N_2);

    for (size_t i = 0; i < N_2; i++)
        result[i] = RandomRange::GetRange(-1, 1);
    return result;
}

void FFT::free_array_float(float* f) {
    pffft_aligned_free(f);
}

float* FFT::fftar_float(PlanMap& plan_map, size_t N, float* in, float* out, float* work, size_t& workN) {
    PFFFT_Setup* plan = read_plan_map_threadsafe(plan_map, N);
    if (N > workN) {
        if (work)
            pffft_aligned_free(work);
        work = (float*)pffft_aligned_malloc(sizeof(float) * N);
        workN = N;
    }
    pffft_transform_ordered(plan, in, out, work, PFFFT_FORWARD);
    out[1] = out[N];
    return work;
}

void FFT::fftsr_destructive_float(PlanMap& plan_map, size_t N, float* in, float* out, float* work) {
    PFFFT_Setup* plan = read_plan_map_threadsafe(plan_map, N);
    pffft_transform_ordered(plan, in, out, work, PFFFT_BACKWARD);
}

void FFT::cleanup(PlanMap& plan_map) {
    auto cleanup_plans = [](PlanMap& plan_map_) {
        for (auto& plan_entry : plan_map_)
            pffft_destroy_setup(plan_entry.second);
        plan_map_.clear();
    };
    cleanup_plans(plan_map);
}
