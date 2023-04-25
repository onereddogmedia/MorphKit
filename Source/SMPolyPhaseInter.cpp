// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMPolyPhaseInter.h"
#include "SMMath.h"

#include <math.h>

using namespace SpectMorph;

using std::vector;

#define OVERSAMPLE 64
#define WIDTH 7

/* set this to at least WIDTH + 2, no problem if it is a little too high */
#define MIN_PADDING 16

#include "SMPolyPhaseCoeffs.h"

static double c_get(int p) {
    if (p >= 0 && p < WIDTH * OVERSAMPLE * 2 - 1)
        return c_7[p];
    else
        return 0;
}

PolyPhaseInter::PolyPhaseInter() {
    /*
     * reorder coefficients, first all coeffs with oversample 0, then oversample 1, ...
     */
    for (int o = 0; o <= OVERSAMPLE; o++) {
        int p = o - 1;

        for (int n = 0; n < WIDTH * 2; n++) {
            x.push_back((float)c_get(p));
            p += OVERSAMPLE;
        }
    }
}

double PolyPhaseInter::get_sample(const Buffer& signal, double pos) {
    const int ipos = (int)pos;

    if (ipos < MIN_PADDING || ipos + MIN_PADDING > int(signal.size())) {
        Buffer shift_signal(MIN_PADDING * 2);

        // shift signal: ipos should be in the center of the generated input signal
        const int shift = (int)shift_signal.size() / 2 - ipos;

        for (int i = 0; i < (int)shift_signal.size(); i++) {
            const int s = i - shift;

            if (s >= 0 && s < (int)signal.size())
                shift_signal[(size_t)i] = signal[(size_t)s];
        }

        return get_sample_no_check(shift_signal, pos + shift);
    } else {
        return get_sample_no_check(signal, pos);
    }
}

double PolyPhaseInter::get_sample_no_check(const Buffer& signal, double pos) {
    const int ipos = (const int)pos;

    const int frac64 = (int)((pos - ipos) * OVERSAMPLE);
    const float frac = (float)((pos - ipos) * OVERSAMPLE - frac64);

    const float* x_a = &x[(size_t)((WIDTH * 2) * (OVERSAMPLE - frac64))];
    const float* x_b = &x[(size_t)((WIDTH * 2) * ((OVERSAMPLE * 2 - frac64 - 1) & (OVERSAMPLE - 1)))];
    const float* s_ptr = &signal[(size_t)(ipos - WIDTH + 1)];

    float result_a = 0, result_b = 0;
    for (int j = 0; j < 2 * WIDTH; j++) {
        result_a += s_ptr[j] * x_a[j];
        result_b += s_ptr[j] * x_b[j];
    }
    return result_a * (1 - frac) + result_b * frac;
}

size_t PolyPhaseInter::get_min_padding() {
    /* Minimum padding before and after the pos sample should be at least MIN_PADDING to use get_sample_no_check() */
    return MIN_PADDING;
}
