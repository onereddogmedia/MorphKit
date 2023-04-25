// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#ifndef SPECTMORPH_MMAP_IN_HH
#define SPECTMORPH_MMAP_IN_HH

#include "SMGenericIn.h"
#include "glib.h"
#include <string>

typedef char GMappedFile;

namespace SpectMorph {

class MMapIn : public GenericIn {
    unsigned char* mapfile;
    unsigned char* mapend;
    unsigned char* pos;

    MMapIn(unsigned char* mapfile, unsigned char* mapend);
    ~MMapIn();

  public:
    static GenericIn* open(const std::string& filename);
    static GenericIn* open_mem(unsigned char* mem_start, unsigned char* mem_end);

    int get_byte(); // like fgetc
    int read(void* ptr, size_t size);
    bool skip(size_t size);
    unsigned char* mmap_mem(size_t& remaining);
    size_t get_pos();
    GenericIn* open_subfile(size_t pos, size_t len);
};

} // namespace SpectMorph

#endif /* SPECTMORPH_MMAP_IN_HH */
