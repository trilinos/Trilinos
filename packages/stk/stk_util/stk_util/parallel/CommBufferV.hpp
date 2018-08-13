// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#ifndef CommBufferV_hpp
#define CommBufferV_hpp

#include <cstddef>
#include <cstring>
#include <vector>

//------------------------------------------------------------------------

namespace stk {

class CommBufferV {
public:
    CommBufferV()
        : data_buffer(), unpack_iterator(data_buffer.begin()) {}
    ~CommBufferV(){}

    size_t size_in_bytes() const { return data_buffer.end() - unpack_iterator; }
    size_t capacity_in_bytes() const { return data_buffer.capacity(); }

    void reserve(size_t num_bytes) {
        data_buffer.reserve(num_bytes);
        unpack_iterator = data_buffer.begin();
    }

    void resize(size_t num_bytes) {
        data_buffer.resize(num_bytes);
        unpack_iterator = data_buffer.begin();
    }

    unsigned char* raw_buffer() { return data_buffer.data(); }
    const unsigned char* raw_buffer() const { return data_buffer.data(); }

    template<typename T>
    void pack(const T& item) {
        pack_internal(&item, 1);
    }

    template<typename T>
    void pack(const T* items, size_t num_items) {
        pack_internal(items, num_items);
    }

    template<typename T>
    void unpack(T& item) {
        enum { item_size_in_bytes = sizeof(T) };
        unsigned char* char_ptr = &(*unpack_iterator);
        T* item_to_unpack = reinterpret_cast<T*>(char_ptr);
        item = *item_to_unpack;
        unpack_iterator += item_size_in_bytes;
    }

    template<typename T>
    void unpack(T* items, size_t num_items) {
        enum { item_size_in_bytes = sizeof(T) };
        unsigned char* char_ptr = &(*unpack_iterator);
        T* items_to_unpack = reinterpret_cast<T*>(char_ptr);
        size_t num_bytes = item_size_in_bytes * num_items;
        std::memcpy(items, items_to_unpack, num_bytes);
        unpack_iterator += item_size_in_bytes * num_items;
    }

private:
    template<typename T>
    void pack_internal(const T* items, size_t num_items) {
        enum { item_size_in_bytes = sizeof(T) };
        const size_t num_bytes = item_size_in_bytes*num_items;
        if (num_bytes > (data_buffer.capacity() - data_buffer.size())) {
            data_buffer.reserve(std::max(num_bytes, data_buffer.capacity()*2));
        }
        const unsigned char* item_chars = reinterpret_cast<const unsigned char*>(items);
        data_buffer.insert(data_buffer.end(), item_chars, item_chars+num_bytes);
        unpack_iterator = data_buffer.begin();
    }

    std::vector<unsigned char> data_buffer;
    std::vector<unsigned char>::iterator unpack_iterator;
};

}

//----------------------------------------------------------------------

#endif

