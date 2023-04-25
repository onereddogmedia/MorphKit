// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMStdioSubIn.h"
#include "SMLeakDebugger.h"

#include <assert.h>
#include <stdio.h>

using namespace SpectMorph;

static LeakDebugger leak_debugger("SpectMorph::StdioSubIn");

GenericIn* StdioSubIn::open(const std::string& filename, size_t pos, size_t len) {
    FILE* file = fopen(filename.c_str(), "rb");

    if (file)
        return new StdioSubIn(file, pos, len);
    else
        return nullptr;
}

StdioSubIn::StdioSubIn(FILE* file_, size_t pos, size_t len) : file(file_) {
    fseek(file, (long)pos, SEEK_SET);
    file_pos = 0;
    file_len = len;

    leak_debugger.add(this);
}

StdioSubIn::~StdioSubIn() {
    assert(file);
    fclose(file);
    leak_debugger.del(this);
}

int StdioSubIn::get_byte() {
    if (file_pos < file_len) {
        int result = fgetc(file);
        file_pos++;
        return result;
    } else {
        return EOF;
    }
}

int StdioSubIn::read(void* ptr, size_t size) {
    if (file_pos < file_len) {
        size_t read_bytes = fread(ptr, 1, size, file);
        file_pos += read_bytes;
        return (int)read_bytes;
    } else {
        return EOF;
    }
}

bool StdioSubIn::skip(size_t size) {
    if (file_pos + size <= file_len) {
        if (fseek(file, (long)size, SEEK_CUR) == 0) {
            file_pos += size;
            return true;
        }
    }
    return false;
}

unsigned char* StdioSubIn::mmap_mem(size_t&) {
    return nullptr;
}

size_t StdioSubIn::get_pos() {
    return file_pos;
}

GenericIn* StdioSubIn::open_subfile(size_t, size_t) {
    assert(false); // implement me
    return nullptr;
}
