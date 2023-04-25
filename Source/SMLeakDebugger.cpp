// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMLeakDebugger.h"
#include "glib.h"
#include <assert.h>

#pragma clang diagnostic ignored "-Wunreachable-code"

using namespace SpectMorph;

using std::map;
using std::string;

#ifdef DEBUG
#    define LEAK (1)
#else
#    define LEAK (0)
#endif

void LeakDebugger::ptr_add(void* p) {
    if (LEAK) {
        std::lock_guard<std::mutex> lock(mutex);

        if (ptr_map[p] != 0)
            printf("LeakDebugger: invalid registration of object type %s detected; ptr_map[p] is %d\n", type.c_str(),
                   ptr_map[p]);

        ptr_map[p]++;
    }
}

void LeakDebugger::ptr_del(void* p) {
    if (LEAK) {
        std::lock_guard<std::mutex> lock(mutex);

        if (ptr_map[p] != 1)
            printf("LeakDebugger: invalid deletion of object type %s detected; ptr_map[p] is %d\n", type.c_str(),
                   ptr_map[p]);

        ptr_map[p]--;
    }
}

LeakDebugger::LeakDebugger(const string& name, std::function<void()> cleanup_function_)
    : type(name), cleanup_function(cleanup_function_) {
}

LeakDebugger::~LeakDebugger() {
    if (LEAK) {
        if (cleanup_function)
            cleanup_function();

        int alive = 0;

        for (map<void*, int>::iterator pi = ptr_map.begin(); pi != ptr_map.end(); pi++) {
            if (pi->second != 0) {
                assert(pi->second == 1);
                alive++;
            }
        }
        if (alive) {
            printf("LeakDebugger (%s) => %d objects remaining\n", type.c_str(), alive);
        }
    }
}
