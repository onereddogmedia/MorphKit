// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_MEM_OUT_HH
#define SPECTMORPH_MEM_OUT_HH

#include "SMGenericOut.h"
#include <string>
#include <vector>

namespace SpectMorph {

class MemOut : public GenericOut {
    std::vector<unsigned char>* output;

  public:
    MemOut(std::vector<unsigned char>* output);
    ~MemOut() override;

    int put_byte(int c) override;
    int write(const void* ptr, size_t size) override;
};

} // namespace SpectMorph

#endif /* SPECTMORPH_MEM_OUT_HH */
