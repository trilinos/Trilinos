// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_STATIC_STORAGE_HPP
#define STOKHOS_STATIC_STORAGE_HPP

#include "Stokhos_StaticArrayTraits.hpp"

#include "Kokkos_Macros.hpp"

#include "Sacado_Traits.hpp"
#include "Stokhos_KokkosTraits.hpp"
#include <sstream>

namespace Stokhos {

  //! Statically allocated storage class
  template <typename ordinal_t, typename value_t, int Num, typename device_t>
  class StaticStorage {
  public:

    static const bool is_static = false;
    static const int static_size = 0;
    static const bool supports_reset = false;

    typedef ordinal_t ordinal_type;
    typedef value_t value_type;
    typedef device_t execution_space;
    typedef value_type& reference;
    typedef volatile value_type& volatile_reference;
    typedef const value_type& const_reference;
    typedef const volatile value_type& const_volatile_reference;
    typedef value_type* pointer;
    typedef volatile value_type* volatile_pointer;
    typedef const value_type* const_pointer;
    typedef const volatile value_type* const_volatile_pointer;
    typedef Stokhos::StaticArrayTraits<value_type,execution_space> ss;

    //! Turn StaticStorage into a meta-function class usable with mpl::apply
    template <typename ord_t, typename val_t = value_t , typename dev_t = device_t >
    struct apply {
      typedef StaticStorage<ord_t,val_t,Num,dev_t> type;
    };

    //! Constructor
    KOKKOS_INLINE_FUNCTION
    StaticStorage(const ordinal_type& sz = 1,
                  const value_type& x = value_type(0.0)) : sz_(sz) {
      ss::fill(coeff_, sz_, x);
    }

    //! Constructor from array
    KOKKOS_INLINE_FUNCTION
    StaticStorage(const ordinal_type& sz, const value_type* x) : sz_(sz) {
      ss::copy(x, coeff_, sz);
    }

    //! Constructor for creating a view (not allowed)
    KOKKOS_INLINE_FUNCTION
    StaticStorage(const ordinal_type& sz, pointer v, bool owned) {}

    //! Copy constructor
    KOKKOS_INLINE_FUNCTION
    StaticStorage(const StaticStorage& s) : sz_(s.sz_) {
      ss::copy(s.coeff_, coeff_, sz_);
    }

    //! Copy constructor
    KOKKOS_INLINE_FUNCTION
    StaticStorage(const volatile StaticStorage& s) : sz_(s.sz_) {
      ss::copy(s.coeff_, coeff_, sz_);
    }

    //! Destructor
    KOKKOS_INLINE_FUNCTION
    ~StaticStorage() {}

    //! Assignment operator
    KOKKOS_INLINE_FUNCTION
    StaticStorage& operator=(const StaticStorage& s) {
      sz_ = s.sz_;
      ss::copy(s.coeff_, coeff_, sz_);
      return *this;
    }

    //! Assignment operator
    KOKKOS_INLINE_FUNCTION
    StaticStorage& operator=(const volatile StaticStorage& s) {
      sz_ = s.sz_;
      ss::copy(s.coeff_, coeff_, sz_);
      return *this;
    }

    //! Assignment operator
    KOKKOS_INLINE_FUNCTION
    volatile StaticStorage& operator=(const StaticStorage& s) volatile {
      sz_ = s.sz_;
      ss::copy(s.coeff_, coeff_, sz_);
      return *this;
    }

    //! Assignment operator
    KOKKOS_INLINE_FUNCTION
    volatile StaticStorage& operator=(const volatile StaticStorage& s) volatile {
      sz_ = s.sz_;
      ss::copy(s.coeff_, coeff_, sz_);
      return *this;
    }

    //! Initialize values to a constant value
    KOKKOS_INLINE_FUNCTION
    void init(const_reference v) {
      ss::fill(coeff_, sz_, v);
    }

    //! Initialize values to a constant value
    KOKKOS_INLINE_FUNCTION
    void init(const_reference v) volatile {
      ss::fill(coeff_, sz_, v);
    }

    //! Initialize values to an array of values
    KOKKOS_INLINE_FUNCTION
    void init(const_pointer v, const ordinal_type& sz = 0) {
      if (sz == 0)
        ss::copy(v, coeff_, sz_);
      else
        ss::copy(v, coeff_, sz);
    }

    //! Initialize values to an array of values
    KOKKOS_INLINE_FUNCTION
    void init(const_pointer v, const ordinal_type& sz = 0) volatile {
      if (sz == 0)
        ss::copy(v, coeff_, sz_);
      else
        ss::copy(v, coeff_, sz);
    }

    //! Load values to an array of values
    KOKKOS_INLINE_FUNCTION
    void load(pointer v) {
      ss::copy(v, coeff_, sz_);
    }

    //! Load values to an array of values
    KOKKOS_INLINE_FUNCTION
    void load(pointer v) volatile {
      ss::copy(v, coeff_, sz_);
    }

    //! Resize to new size (values are preserved)
    KOKKOS_INLINE_FUNCTION
    void resize(const ordinal_type& sz) {
      if (sz > sz_)
        ss::fill(coeff_+sz_, sz-sz_, value_type(0.0));
      sz_ = sz;
    }

    //! Resize to new size (values are preserved)
    KOKKOS_INLINE_FUNCTION
    void resize(const ordinal_type& sz) volatile {
      if (sz > sz_)
        ss::fill(coeff_+sz_, sz-sz_, value_type(0.0));
      sz_ = sz;
    }

    //! Reset storage to given array, size, and stride
    KOKKOS_INLINE_FUNCTION
    void shallowReset(pointer v, const ordinal_type& sz,
                      const ordinal_type& stride, bool owned) {}

    //! Reset storage to given array, size, and stride
    KOKKOS_INLINE_FUNCTION
    void shallowReset(pointer v, const ordinal_type& sz,
                      const ordinal_type& stride, bool owned) volatile {}

    //! Return size
    KOKKOS_INLINE_FUNCTION
    ordinal_type size() const { return sz_; }

     //! Return size
    KOKKOS_INLINE_FUNCTION
    ordinal_type size() const volatile { return sz_; }

    //! Coefficient access (avoid if possible)
    KOKKOS_INLINE_FUNCTION
    const_reference operator[] (const ordinal_type& i) const {
      return coeff_[i];
    }

    //! Coefficient access (avoid if possible)
    KOKKOS_INLINE_FUNCTION
    const_volatile_reference operator[] (const ordinal_type& i) const volatile {
      return coeff_[i];
    }

    //! Coefficient access (avoid if possible)
    KOKKOS_INLINE_FUNCTION
    reference operator[] (const ordinal_type& i) { return coeff_[i]; }

    //! Coefficient access (avoid if possible)
    KOKKOS_INLINE_FUNCTION
    volatile_reference operator[] (const ordinal_type& i) volatile {
      return coeff_[i];
    }

     template <int i>
    KOKKOS_INLINE_FUNCTION
    reference getCoeff() { return coeff_[i]; }

    template <int i>
    KOKKOS_INLINE_FUNCTION
    volatile_reference getCoeff() volatile { return coeff_[i]; }

    template <int i>
    KOKKOS_INLINE_FUNCTION
    const_volatile_reference getCoeff() const volatile { return coeff_[i]; }

     template <int i>
    KOKKOS_INLINE_FUNCTION
    const_reference getCoeff() const { return coeff_[i]; }

    //! Get coefficients
    KOKKOS_INLINE_FUNCTION
    const_volatile_pointer coeff() const volatile { return coeff_; }

    //! Get coefficients
    KOKKOS_INLINE_FUNCTION
    const_pointer coeff() const { return coeff_; }

    //! Get coefficients
    KOKKOS_INLINE_FUNCTION
    volatile_pointer coeff() volatile { return coeff_; }

    //! Get coefficients
    KOKKOS_INLINE_FUNCTION
    pointer coeff() { return coeff_; }

  private:

    //! Coefficient values
    value_type coeff_[Num];

    //! Size of array used
    ordinal_type sz_;

  };

}

#include "Stokhos_StorageHelpers.hpp"
STOKHOS_STORAGE_HELPER_STRINGNAME_STATIC(StaticStorage)

#endif // STOKHOS_STATIC_STORAGE_HPP
