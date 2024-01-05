// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMMemOut.h"
#include "SMLeakDebugger.h"
#include <stdio.h>

using namespace SpectMorph;

using std::vector;

static LeakDebugger leak_debugger("SpectMorph::MemOut");

MemOut::MemOut(vector<unsigned char>* output_) : output(output_) {
    leak_debugger.add(this);
}

MemOut::~MemOut() {
    leak_debugger.del(this);
}

int MemOut::put_byte(int c) {
    output->push_back((unsigned char)c);
    return c;
}

int MemOut::write(const void* ptr, size_t size) {
    const unsigned char* uptr = reinterpret_cast<const unsigned char*>(ptr);
    output->insert(output->end(), uptr, uptr + size);
    return (int)size;
}
