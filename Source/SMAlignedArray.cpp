// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMAlignedArray.h"

#include "glib.h"

namespace SpectMorph {

/* --- memory utils --- */
void* malloc_aligned(size_t total_size, size_t alignment, uint8** free_pointer) {
    const bool alignment_power_of_2 = (alignment & (alignment - 1)) == 0;
    const size_t cache_line_size = 64; // ensure that no false sharing will occur (at begin and end of data)
    if (alignment_power_of_2) {
        // for power of 2 alignment, we guarantee also cache line alignment
        alignment = std::max(alignment, cache_line_size);
        uint8* aligned_mem = (uint8*)malloc(total_size + (alignment - 1) + (cache_line_size - 1));
        *free_pointer = aligned_mem;
        if ((std::ptrdiff_t)aligned_mem % (std::ptrdiff_t)alignment)
            aligned_mem += (std::ptrdiff_t)alignment - (std::ptrdiff_t)aligned_mem % (std::ptrdiff_t)alignment;
        return aligned_mem;
    } else {
        uint8* aligned_mem = (uint8*)malloc(total_size + (alignment - 1) + (cache_line_size - 1) * 2);
        *free_pointer = aligned_mem;
        if ((std::ptrdiff_t)aligned_mem % (std::ptrdiff_t)cache_line_size)
            aligned_mem +=
                (std::ptrdiff_t)cache_line_size - (std::ptrdiff_t)aligned_mem % (std::ptrdiff_t)cache_line_size;
        if ((std::ptrdiff_t)aligned_mem % (std::ptrdiff_t)alignment)
            aligned_mem += (std::ptrdiff_t)alignment - (std::ptrdiff_t)aligned_mem % (std::ptrdiff_t)alignment;
        return aligned_mem;
    }
}

} // namespace SpectMorph
