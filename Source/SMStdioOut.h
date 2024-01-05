// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_STDIO_OUT_HH
#define SPECTMORPH_STDIO_OUT_HH

#include "SMGenericOut.h"
#include <string>

namespace SpectMorph {

class StdioOut : public GenericOut {
    FILE* file;

    StdioOut(FILE* file);

  public:
    static GenericOut* open(const std::string& filename);

    ~StdioOut() override;
    int put_byte(int c) override;
    int write(const void* ptr, size_t size) override;
};

} // namespace SpectMorph

#endif /* SPECTMORPH_STDIO_OUT_HH */
