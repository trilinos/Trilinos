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

#ifndef PACKAGES_STK_STK_UTIL_STK_UTIL_UTIL_STKVECTOR_HPP_
#define PACKAGES_STK_STK_UTIL_STK_UTIL_UTIL_STKVECTOR_HPP_

#include <Kokkos_Core.hpp>
#include <stk_util/stk_kokkos_macros.h>

namespace stk
{

template <typename Datatype>
class Vector
{
public:
    Vector(const std::string &n) : mSize(0), deviceVals(n), hostVals(Kokkos::create_mirror_view(deviceVals))
    {
    }
    Vector() : Vector(get_default_name())
    {
    }
    Vector(const std::string &n, size_t s) : mSize(s), deviceVals(n, mSize), hostVals(Kokkos::create_mirror_view(deviceVals))
    {
    }
    Vector(size_t s) : Vector(get_default_name(), s)
    {
    }
    Vector(const std::string &n, size_t s, Datatype init) : Vector(n, s)
    {
        Kokkos::deep_copy(hostVals, init);
    }
    Vector(size_t s, Datatype init) : Vector(get_default_name(), s, init)
    {
    }
    ~Vector() {}

    std::string name() const { return hostVals.label(); }

    STK_FUNCTION size_t size() const { return mSize; }
    STK_FUNCTION bool empty() const { return mSize == 0; }
    size_t capacity() const
    {
        return hostVals.size();
    }

    void resize(size_t s)
    {
        resize(s, 0);
    }
    void resize(size_t s, Datatype init)
    {
        if(s > capacity())
            grow_to_size(s);
        if(s > mSize)
            for(size_t i=mSize; i<s; i++)
                hostVals(i) = init;
        mSize = s;
    }

    void clear()
    {
        mSize = 0;
    }

    Datatype & operator[](size_t i) const
    {
        return hostVals(i);
    }
    STK_FUNCTION Datatype & device_get(size_t i) const
    {
        return deviceVals(i);
    }

    void push_back(Datatype val)
    {
        if(mSize >= capacity())
            grow_to_size(mSize + get_push_back_increase_size());
        hostVals[mSize] = val;
        mSize++;
    }

    void copy_host_to_device()
    {
        Kokkos::deep_copy(deviceVals, hostVals);
    }
    void copy_device_to_host()
    {
        Kokkos::deep_copy(hostVals, deviceVals);
    }

protected:
#ifdef KOKKOS_ENABLE_CUDA
  using DeviceSpace = Kokkos::CudaSpace;
#else
  using DeviceSpace = Kokkos::HostSpace;
#endif
  typedef Kokkos::View<Datatype *, DeviceSpace> DeviceType;
    typedef typename DeviceType::HostMirror HostType;

    virtual DeviceType get_new_vals_of_size(size_t s)
    {
        return DeviceType(deviceVals.label(), s);
    }

private:
    size_t mSize;
    DeviceType deviceVals;
    HostType hostVals;

    const char *get_default_name() { return "UnnamedStkVector"; }

    size_t get_push_back_increase_size() const
    {
        if(mSize == 0)
            return 1;
        return mSize;
    }

    void grow_to_size(size_t s)
    {
        deviceVals = get_new_vals_of_size(s);
        HostType tmp = Kokkos::create_mirror_view(deviceVals);
        copy_into_bigger(tmp, hostVals);
        hostVals = tmp;
    }

    void copy_into_bigger(HostType& dst, HostType& src)
    {
        for(size_t i = 0; i < src.size(); i++)
            dst[i] = src[i];
    }
};

}

#endif /* PACKAGES_STK_STK_UTIL_STK_UTIL_UTIL_STKVECTOR_HPP_ */
