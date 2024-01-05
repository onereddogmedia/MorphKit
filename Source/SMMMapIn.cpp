// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMMMapIn.h"
#include "SMLeakDebugger.h"
#include "SMUtils.h"

#include <assert.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

using namespace SpectMorph;

static LeakDebugger leak_debugger("SpectMorph::MMapIn");

GenericIn* MMapIn::open(const std::string&) {
    return nullptr;
}

MMapIn::MMapIn(const unsigned char* mapfile_, const unsigned char* mapend_) : mapfile(mapfile_), mapend(mapend_) {
    pos = static_cast<const unsigned char*>(mapfile);

    leak_debugger.add(this);
}

MMapIn::~MMapIn() {
    leak_debugger.del(this);
}

int MMapIn::get_byte() {
    if (pos < mapend)
        return *pos++;
    else
        return EOF;
}

int MMapIn::read(void* ptr, size_t size) {
    if (pos + size <= mapend) {
        memcpy(ptr, pos, size);
        pos += size;
        return (int)size;
    } else
        return 0;
}

bool MMapIn::skip(size_t size) {
    if (pos + size <= mapend) {
        pos += size;
        return true;
    } else {
        return false;
    }
}

const unsigned char* MMapIn::mmap_mem(size_t& remaining) {
    remaining = (size_t)(mapend - pos);
    return pos;
}

size_t MMapIn::get_pos() {
    return (size_t)(pos - mapfile);
}

GenericIn* MMapIn::open_subfile(size_t pos_, size_t len) {
    return new MMapIn(mapfile + pos_, mapfile + pos_ + len);
}
