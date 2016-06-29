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
    Vector(const std::string &n, size_t s, Datatype init) : mSize(s), deviceVals(n, mSize), hostVals(Kokkos::create_mirror_view(deviceVals))
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
    typedef Kokkos::View<Datatype*> DeviceType;
    typedef typename DeviceType::HostMirror HostType;

    virtual HostType get_new_vals_of_size(size_t s)
    {
        return HostType(hostVals.label(), s);
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
        HostType tmp = get_new_vals_of_size(s);
        copy_into_bigger(tmp, hostVals);
        hostVals = tmp;
        hostVals = Kokkos::create_mirror_view(hostVals);
    }

    void copy_into_bigger(HostType& dst, HostType& src)
    {
        for(size_t i = 0; i < src.size(); i++)
            dst[i] = src[i];
    }
};

}

#endif /* PACKAGES_STK_STK_UTIL_STK_UTIL_UTIL_STKVECTOR_HPP_ */
