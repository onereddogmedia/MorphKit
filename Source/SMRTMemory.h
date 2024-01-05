// Licensed GNU LGPL v2.1 or later: http://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include <cassert>
#include <vector>

#include "SMAudio.h"
#include "SMMath.h"
#include "SMUtils.h"

namespace SpectMorph {

class RTMemoryArea {
    SPECTMORPH_CLASS_NON_COPYABLE(RTMemoryArea);

    std::vector<unsigned char> m_mem;
    std::vector<void*> m_malloc_mem;
    size_t m_used = 0;
    static constexpr size_t initial_size = 1024 * 1024;
    static constexpr size_t grow_size = 32 * 1024;

  public:
    RTMemoryArea() : m_mem(initial_size) {
    }
    void* alloc(size_t n) {
        n = ((n + 63) / 64) * 64;

        void* result;
        if (m_used + n > m_mem.size()) {
            /* if we don't have enough memory available, fall back to malloc */
            result = malloc(n);
            m_malloc_mem.push_back(result);
        } else {
            result = m_mem.data() + m_used;
        }
        m_used += n;
        return result;
    }
    void free_all() {
        for (auto ptr : m_malloc_mem)
            free(ptr);
        m_malloc_mem.clear();

        if (m_used > m_mem.size()) {
            /* if we didn't have enough memory last time, make the memory area larger
             *
             * this really should *never* happen (the initial size should be large
             * enough), but we do the fallback-to-malloc to not crash in cases
             * where the initial size was too small
             */
            m_mem.resize(m_used + grow_size);
        }

        m_used = 0;
    }
};

template <class T> class RTVector {
    SPECTMORPH_CLASS_NON_COPYABLE(RTVector);

    RTMemoryArea* m_memory_area = nullptr;
    T* m_start = nullptr;
    size_t m_size = 0;
    size_t m_capacity = 0;

  public:
    RTVector(RTMemoryArea* memory_area) : m_memory_area(memory_area) {
    }
    void assign(const std::vector<T>& vec) {
        assert(m_size == 0 && m_capacity == 0);

        set_capacity(vec.size());
        std::copy(vec.begin(), vec.end(), m_start);
        m_size = vec.size();
    }
    void assign(const RTVector<T>& vec) {
        assert(m_size == 0 && m_capacity == 0);

        set_capacity(vec.size());
        std::copy(vec.m_start, vec.m_start + vec.m_size, m_start);
        m_size = vec.size();
    }
    size_t size() const {
        return m_size;
    }
    void set_capacity(size_t capacity) {
        assert(m_size == 0 && m_capacity == 0);

        m_start = (T*)m_memory_area->alloc(sizeof(T) * capacity);
        m_capacity = capacity;
    }
    void push_back(const T& t) {
        assert(m_size < m_capacity);
        m_start[m_size++] = t;
    }
    T& back() {
        return m_start[m_size - 1];
    }
    T& operator[](size_t idx) {
        return m_start[idx];
    }
    const T& operator[](size_t idx) const {
        return m_start[idx];
    }
};

class RTAudioBlock {
  public:
    RTAudioBlock(RTMemoryArea* memory_area) : freqs(memory_area), mags(memory_area), noise(memory_area) {
    }
    void assign(const RTAudioBlock& audio_block) {
        freqs.assign(audio_block.freqs);
        mags.assign(audio_block.mags);
        noise.assign(audio_block.noise);
    }
    void assign(const AudioBlock& audio_block) {
        freqs.assign(audio_block.freqs);
        mags.assign(audio_block.mags);
        noise.assign(audio_block.noise);
    }
    RTVector<uint16_t> freqs;
    RTVector<uint16_t> mags;
    RTVector<uint16_t> noise;

    double freqs_f(size_t i) const {
        return sm_ifreq2freq(freqs[i]);
    }

    double mags_f(size_t i) const {
        return sm_idb2factor(mags[i]);
    }

    double noise_f(size_t i) const {
        return sm_idb2factor(noise[i]);
    }

    void sort_freqs();
};

} // namespace SpectMorph
