// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
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
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
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

#include "Kokkos_Core.hpp"

namespace stk
{

template <typename Datatype>
class NgpVector
{
  using HostSpace = Kokkos::DefaultHostExecutionSpace;
public:
    NgpVector(const std::string &n) : NgpVector(n, 0)
    {
    }
    NgpVector() : NgpVector(get_default_name())
    {
    }
    NgpVector(const std::string &n, size_t s)
        : mSize(s),
          deviceVals(Kokkos::view_alloc(Kokkos::WithoutInitializing, n), mSize),
          hostVals(Kokkos::create_mirror_view(Kokkos::WithoutInitializing, deviceVals))
        {
    }
    NgpVector(size_t s) : NgpVector(get_default_name(), s)
    {
    }
    NgpVector(const std::string &n, size_t s, Datatype init) : NgpVector(n, s)
    {
        Kokkos::deep_copy(hostVals, init);
    }
    NgpVector(size_t s, Datatype init) : NgpVector(get_default_name(), s, init) {}

    std::string name() const { return hostVals.label(); }

    auto view_host() { return Kokkos::subview(hostVals, Kokkos::make_pair(size_t(0), size())); }
    auto view_device() { return Kokkos::subview(deviceVals, Kokkos::make_pair(size_t(0), size())); }

    KOKKOS_FUNCTION size_t size() const { return mSize; }
    KOKKOS_FUNCTION bool empty() const { return mSize == 0; }
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
    KOKKOS_FUNCTION Datatype & device_get(size_t i) const
    {
        return deviceVals(i);
    }

protected:
#ifdef KOKKOS_ENABLE_CUDA
  using DeviceSpace = Kokkos::CudaSpace;
#elif defined(KOKKOS_ENABLE_HIP)
  using DeviceSpace = Kokkos::HIPSpace;
#else
  using DeviceSpace = Kokkos::HostSpace;
#endif
public:
    template <class Device>
    KOKKOS_FUNCTION Datatype & get(
      typename std::enable_if<
        std::is_same<typename Device::execution_space, DeviceSpace::execution_space>::value,
        size_t>::type i) const
    {
      return deviceVals(i);
    }
#ifdef STK_ENABLE_GPU
    template <class Device>
    KOKKOS_FUNCTION Datatype & get(
      typename std::enable_if<
        !std::is_same<typename Device::execution_space, DeviceSpace::execution_space>::value,
        size_t>::type i) const
    {
      return hostVals(i);
    }
#endif

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

    static const char *get_default_name() { return "UnnamedStkVector"; }

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
