#pragma once

#include <assert.h>

#define FALSE 0
#define TRUE 1

#define g_return_if_fail(x)                                                                                            \
    if (!(x)) {                                                                                                        \
        return;                                                                                                        \
    }

#define g_return_val_if_fail(x, y)                                                                                     \
    if (!(x)) {                                                                                                        \
        return (y);                                                                                                    \
    }

#define g_assert_not_reached() assert(false)
#define g_assert(x) assert(x)

#define g_warning printf
#define g_printerr printf
