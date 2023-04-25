// Licensed GNU LGPL v3 or later: http://www.gnu.org/licenses/lgpl.html

#include "SMStdioIn.h"
#include "SMLeakDebugger.h"
#include "SMStdioSubIn.h"

#include <assert.h>
#include <stdio.h>

using namespace SpectMorph;

using std::string;

static LeakDebugger leak_debugger("SpectMorph::StdioIn");

GenericIn* StdioIn::open(const string& filename) {
    FILE* file = fopen(filename.c_str(), "rb");

    if (file) {
        return new StdioIn(file, filename);
    } else {
        printf("failed to open: %s\n", filename.c_str());
        return nullptr;
    }
}

StdioIn::StdioIn(FILE* f, const string& filename_) : file(f), filename(filename_) {
    leak_debugger.add(this);
}

StdioIn::~StdioIn() {
    assert(file);
    fclose(file);
    leak_debugger.del(this);
}

int StdioIn::get_byte() {
    return fgetc(file);
}

int StdioIn::read(void* ptr, size_t size) {
    return (int)fread(ptr, 1, size, file);
}

bool StdioIn::skip(size_t size) {
    if (fseek(file, (long)size, SEEK_CUR) == 0)
        return true;
    return false;
}

unsigned char* StdioIn::mmap_mem(size_t&) {
    return nullptr;
}

size_t StdioIn::get_pos() {
    return (size_t)ftell(file);
}

GenericIn* StdioIn::open_subfile(size_t pos, size_t len) {
    return StdioSubIn::open(filename, pos, len);
}
