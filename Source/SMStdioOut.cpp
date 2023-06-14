// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMStdioOut.h"
#include "SMLeakDebugger.h"
#include <stdio.h>

using namespace SpectMorph;

static LeakDebugger leak_debugger("SpectMorph::StdioOut");

GenericOut* StdioOut::open(const std::string& filename) {
    FILE* file = fopen(filename.c_str(), "wb");

    if (file)
        return new StdioOut(file);
    else
        return nullptr;
}

StdioOut::StdioOut(FILE* file_) : file(file_) {
    leak_debugger.add(this);
}

StdioOut::~StdioOut() {
    if (file != nullptr) {
        fclose(file);
        file = nullptr;
    }
    leak_debugger.del(this);
}

int StdioOut::put_byte(int c) {
    return fputc(c, file);
}

size_t StdioOut::write(const void* ptr, size_t size) {
    if (size == 0)
        return 0;
    return fwrite(ptr, 1, size, file);
}
