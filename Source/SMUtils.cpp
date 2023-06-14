// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMUtils.h"

#include <codecvt>
#include <locale>
#include <string>

#include "glib.h"
#include <stdarg.h>
#include <sys/stat.h>

#ifdef SM_OS_MACOS
#    include <sys/time.h>
#    include <xlocale.h>
#endif

#ifdef SM_OS_LINUX
#    include <locale.h>
#    include <sys/time.h>
#endif

#ifdef SM_OS_WINDOWS
#    include "Windows.h"
#    include "shlobj.h"
#    include <sys/timeb.h>

int gettimeofday(struct timeval* tp, struct timezone*) {
    // Note: some broken versions only have 8 trailing zero's, the correct epoch has 9 trailing zero's
    static const uint64_t EPOCH = ((uint64_t)116444736000000000ULL);

    SYSTEMTIME system_time;
    FILETIME file_time;
    uint64_t time;

    GetSystemTime(&system_time);
    SystemTimeToFileTime(&system_time, &file_time);
    time = ((uint64_t)file_time.dwLowDateTime);
    time += ((uint64_t)file_time.dwHighDateTime) << 32;

    tp->tv_sec = (long)((time - EPOCH) / 10000000L);
    tp->tv_usec = (long)(system_time.wMilliseconds * 1000);
    return 0;
}

int vasprintf(char** strp, const char* fmt, va_list ap) {
    va_list ap2;
    va_copy(ap2, ap);
    char tmp[1];
    int size = vsnprintf(tmp, 1, fmt, ap2);
    if (size <= 0) {
        va_end(ap2);
        return size;
    }
    va_end(ap2);
    size += 1;
    *strp = (char*)malloc(size * sizeof(char));
    return vsnprintf(*strp, size, fmt, ap);
}

#endif

using std::string;
using std::u32string;
using std::vector;

static string string_current_vprintf(const char* format, va_list vargs) {
    string s;
    char* str = nullptr;
    if (vasprintf(&str, format, vargs) >= 0 && str) {
        s = str;
        free(str);
    } else {
        s = format;
        free(str);
    }
    return s;
}

namespace SpectMorph {

string string_vprintf(const char* format, va_list vargs) {
    return string_current_vprintf(format, vargs);
}

string string_printf(const char* format, ...) {
    string str;
    va_list args;
    va_start(args, format);
    str = string_vprintf(format, args);
    va_end(args);
    return str;
}

string string_locale_printf(const char* format, ...) {
    string str;
    va_list args;
    va_start(args, format);
    str = string_current_vprintf(format, args);
    va_end(args);
    return str;
}

void sm_printf(const char* format, ...) {
    string str;
    va_list args;
    va_start(args, format);
    str = string_vprintf(format, args);
    va_end(args);

    printf("%s", str.c_str());
}

static string pkg_data_dir = "";
static string icloud_data_dir = "";

void sm_set_pkg_data_dir(const string& data_dir) {
    pkg_data_dir = data_dir;
}

void sm_set_icloud_data_dir(const string& data_dir) {
    icloud_data_dir = data_dir;
}

std::string sm_get_install_dir(InstallDir p) {
    switch (p) {
        case INSTALL_DIR_INSTRUMENTS:
            return pkg_data_dir;
        case INSTALL_DIR_ICLOUD:
            return icloud_data_dir;
    }
    return "";
}

string sha1_hash(const unsigned char*, size_t) {
    //  char *result = g_compute_checksum_for_data (G_CHECKSUM_SHA1, data, len);
    //  string hash = result;
    //  free (result);
    string hash = "morphwizpro";

    return hash;
}

string sha1_hash(const string& str) {
    return sha1_hash(reinterpret_cast<const unsigned char*>(str.data()), str.size());
}

double get_time() {
    /* return timestamp in seconds as double */
    timeval tv;
    gettimeofday(&tv, nullptr);

    return double(tv.tv_sec) + double(tv.tv_usec) / 1000000.0;
}

#ifdef __clang__
#    pragma clang diagnostic push
#    pragma clang diagnostic ignored "-Wdeprecated-declarations"
#endif

string to_utf8(const u32string& str) {
    std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> conv;
    return conv.to_bytes(str);
}

u32string to_utf32(const string& utf8) {
    std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> conv;
    return conv.from_bytes(utf8);
}

#ifdef __clang__
#    pragma clang diagnostic pop
#endif

double sm_atof(const char* str) // always use . as decimal seperator
{
    return atof(str);
}

} // namespace SpectMorph
