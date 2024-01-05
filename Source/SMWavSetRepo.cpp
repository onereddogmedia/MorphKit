// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#include "SMWavSetRepo.h"

using namespace SpectMorph;

using std::string;

WavSetRepo::WavSetRepo() {
    wav_set_map.reset(new lru_cache_using_std<std::string, std::shared_ptr<WavSet>, std::map>(open, 4));
}

WavSetRepo::~WavSetRepo() {
}

WavSet* WavSetRepo::get(const string& filename) {
    std::lock_guard<std::mutex> lock(mutex);
    auto wav = (*wav_set_map)(filename);
    return wav.get();
}

void WavSetRepo::reload() {
    wav_set_map->reload();
}

std::shared_ptr<WavSet> WavSetRepo::open(const std::string& s) {
    std::shared_ptr<WavSet> wav_set(new WavSet());
    Error err = wav_set->load(s);
    if (err.code() != Error::Code::NONE) {
        printf("failed to read: %s\n", s.c_str());
        wav_set = nullptr;
    }
    return wav_set;
}
