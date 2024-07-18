// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_DYNAMIC_STRIDED_STORAGE_HPP
#define STOKHOS_DYNAMIC_STRIDED_STORAGE_HPP

#include "Stokhos_DynArrayTraits.hpp"

#include "Kokkos_Macros.hpp"

#include "Sacado_Traits.hpp"
#include "Stokhos_KokkosTraits.hpp"
#include <sstream>

namespace Stokhos {

  template <typename ordinal_t, typename value_t, typename device_t>
  class DynamicStridedStorage {
  public:

    static const bool is_static = false;
    static const int static_size = 0;
    static const bool supports_reset = true;

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
    typedef Stokhos::DynArrayTraits<value_type,execution_space> ds;

    //! Turn DynamicStridedStorage into a meta-function class usable with mpl::apply
    template <typename ord_t, typename val_t = value_t, typename dev_t = device_t >
    struct apply {
      typedef DynamicStridedStorage<ord_t,val_t,dev_t> type;
    };

    //! Constructor
    KOKKOS_INLINE_FUNCTION
    DynamicStridedStorage(const ordinal_type& sz = 1,
                          const value_type& x = value_type(0.0)) :
      sz_(sz), stride_(1), is_owned_(true) {
      coeff_ = ds::get_and_fill(sz_, x);
    }

    //! Constructor from array
    KOKKOS_INLINE_FUNCTION
    DynamicStridedStorage(const ordinal_type& sz, const value_type* x) :
      sz_(sz), stride_(1), is_owned_(true) {
      coeff_ = ds::get_and_fill(x, sz_);
    }

    //! Constructor for creating a view
    KOKKOS_INLINE_FUNCTION
    DynamicStridedStorage(const ordinal_type& sz, pointer v, bool owned) :
      coeff_(v), sz_(sz), stride_(1), is_owned_(owned) {}

    //! Constructor
    KOKKOS_INLINE_FUNCTION
    DynamicStridedStorage(const DynamicStridedStorage& s) :
    sz_(s.sz_), stride_(1), is_owned_(true) {
      if (s.stride_ == 1)
        coeff_ = ds::get_and_fill(s.coeff_, sz_);
      else {
        coeff_ = ds::get_and_fill(sz_);
        for (ordinal_type i=0; i<sz_; ++i)
          coeff_[i] = s[i];
      }
    }

    //! Constructor
    KOKKOS_INLINE_FUNCTION
    DynamicStridedStorage(const volatile DynamicStridedStorage& s) :
    sz_(s.sz_), stride_(1), is_owned_(true) {
      if (s.stride_ == 1)
        coeff_ = ds::get_and_fill(s.coeff_, sz_);
      else {
        coeff_ = ds::get_and_fill(sz_);
        for (ordinal_type i=0; i<sz_; ++i)
          coeff_[i] = s[i];
      }
    }

    //! Destructor
    KOKKOS_INLINE_FUNCTION
    ~DynamicStridedStorage() {
      if (is_owned_) ds::destroy_and_release(coeff_, sz_*stride_);
    }

    //! Assignment operator
    KOKKOS_INLINE_FUNCTION
    DynamicStridedStorage& operator=(const DynamicStridedStorage& s) {
      if (&s != this) {
        if (s.sz_ != sz_) {
          if (is_owned_)
            ds::destroy_and_release(coeff_, sz_*stride_);
          if (s.stride_ == 1)
            coeff_ = ds::get_and_fill(s.coeff_, s.sz_);
          else {
            coeff_ = ds::get_and_fill(s.sz_);
            for (ordinal_type i=0; i<s.sz_; ++i)
              coeff_[i] = s[i];
          }
          sz_ = s.sz_;
          stride_ = 1;
          is_owned_ = true;
        }
        else {
          if (stride_ == 1 and s.stride_ == 1)
            ds::copy(s.coeff_, coeff_, sz_);
          else
            for (ordinal_type i=0; i<s.sz_; ++i)
              coeff_[i*stride_] = s[i];
        }
      }
      return *this;
    }

    //! Assignment operator
    KOKKOS_INLINE_FUNCTION
    DynamicStridedStorage& operator=(const volatile DynamicStridedStorage& s) {
      if (&s != this) {
        if (s.sz_ != sz_) {
          if (is_owned_)
            ds::destroy_and_release(coeff_, sz_*stride_);
          if (s.stride_ == 1)
            coeff_ = ds::get_and_fill(s.coeff_, s.sz_);
          else {
            coeff_ = ds::get_and_fill(s.sz_);
            for (ordinal_type i=0; i<s.sz_; ++i)
              coeff_[i] = s[i];
          }
          sz_ = s.sz_;
          stride_ = 1;
          is_owned_ = true;
        }
        else {
          if (stride_ == 1 and s.stride_ == 1)
            ds::copy(s.coeff_, coeff_, sz_);
          else
            for (ordinal_type i=0; i<s.sz_; ++i)
              coeff_[i*stride_] = s[i];
        }
      }
      return *this;
    }

    //! Assignment operator
    KOKKOS_INLINE_FUNCTION
    volatile DynamicStridedStorage&
    operator=(const DynamicStridedStorage& s) volatile {
      if (&s != this) {
        if (s.sz_ != sz_) {
          if (is_owned_)
            ds::destroy_and_release(coeff_, sz_*stride_);
          if (s.stride_ == 1)
            coeff_ = ds::get_and_fill(s.coeff_, s.sz_);
          else {
            coeff_ = ds::get_and_fill(s.sz_);
            for (ordinal_type i=0; i<s.sz_; ++i)
              coeff_[i] = s[i];
          }
          sz_ = s.sz_;
          stride_ = 1;
          is_owned_ = true;
        }
        else {
          if (stride_ == 1 and s.stride_ == 1)
            ds::copy(s.coeff_, coeff_, sz_);
          else
            for (ordinal_type i=0; i<s.sz_; ++i)
              coeff_[i*stride_] = s[i];
        }
      }
      return *this;
    }

    //! Assignment operator
    KOKKOS_INLINE_FUNCTION
    volatile DynamicStridedStorage&
    operator=(const volatile DynamicStridedStorage& s) volatile {
      if (&s != this) {
        if (s.sz_ != sz_) {
          if (is_owned_)
            ds::destroy_and_release(coeff_, sz_*stride_);
          if (s.stride_ == 1)
            coeff_ = ds::get_and_fill(s.coeff_, s.sz_);
          else {
            coeff_ = ds::get_and_fill(s.sz_);
            for (ordinal_type i=0; i<s.sz_; ++i)
              coeff_[i] = s[i];
          }
          sz_ = s.sz_;
          stride_ = 1;
          is_owned_ = true;
        }
        else {
          if (stride_ == 1 and s.stride_ == 1)
            ds::copy(s.coeff_, coeff_, sz_);
          else
            for (ordinal_type i=0; i<s.sz_; ++i)
              coeff_[i*stride_] = s[i];
        }
      }
      return *this;
    }

    //! Initialize values to a constant value
    KOKKOS_INLINE_FUNCTION
    void init(const_reference v) {
      if (stride_ == 1)
        ds::fill(coeff_, sz_, v);
      else
        for (ordinal_type i=0; i<sz_; ++i)
          coeff_[i*stride_] = v;
    }

    //! Initialize values to a constant value
    KOKKOS_INLINE_FUNCTION
    void init(const_reference v) volatile {
      if (stride_ == 1)
        ds::fill(coeff_, sz_, v);
      else
        for (ordinal_type i=0; i<sz_; ++i)
          coeff_[i*stride_] = v;
    }

    //! Initialize values to an array of values
    KOKKOS_INLINE_FUNCTION
    void init(const_pointer v, const ordinal_type& sz = 0) {
      ordinal_type my_sz = sz;
      if (sz == 0)
        my_sz = sz_;
      if (stride_ == 1)
        ds::copy(v, coeff_, my_sz);
      else
        for (ordinal_type i=0; i<my_sz; ++i)
          coeff_[i*stride_] = v[i];
    }

    //! Initialize values to an array of values
    KOKKOS_INLINE_FUNCTION
    void init(const_pointer v, const ordinal_type& sz = 0) volatile {
      ordinal_type my_sz = sz;
      if (sz == 0)
        my_sz = sz_;
      if (stride_ == 1)
        ds::copy(v, coeff_, my_sz);
      else
        for (ordinal_type i=0; i<my_sz; ++i)
          coeff_[i*stride_] = v[i];
    }

    //! Load values to an array of values
    KOKKOS_INLINE_FUNCTION
    void load(pointer v) {
      if (stride_ == 1)
        copy(v, coeff_, sz_);
      for (ordinal_type i=0; i<sz_; ++i)
        coeff_[i*stride_] = v[i];
    }

    //! Load values to an array of values
    KOKKOS_INLINE_FUNCTION
    void load(pointer v) volatile {
      if (stride_ == 1)
        copy(v, coeff_, sz_);
      for (ordinal_type i=0; i<sz_; ++i)
        coeff_[i*stride_] = v[i];
    }

    //! Resize to new size (values are preserved)
    KOKKOS_INLINE_FUNCTION
    void resize(const ordinal_type& sz) {
      if (sz != sz_) {
        value_type *coeff_new = ds::get_and_fill(sz);
        ordinal_type my_sz = sz_;
        if (sz_ > sz)
          my_sz = sz;
        if (stride_ == 1)
          ds::copy(coeff_, coeff_new, my_sz);
        else
          for (ordinal_type i=0; i<my_sz; ++i)
            coeff_new[i] = coeff_[i*stride_];
        if (is_owned_)
          ds::destroy_and_release(coeff_, sz_*stride_);
        coeff_ = coeff_new;
        sz_ = sz;
        stride_ = 1;
        is_owned_ = true;
      }
    }

    //! Resize to new size (values are preserved)
    KOKKOS_INLINE_FUNCTION
    void resize(const ordinal_type& sz) volatile {
      if (sz != sz_) {
        value_type *coeff_new = ds::get_and_fill(sz);
        ordinal_type my_sz = sz_;
        if (sz_ > sz)
          my_sz = sz;
        if (stride_ == 1)
          ds::copy(coeff_, coeff_new, my_sz);
        else
          for (ordinal_type i=0; i<my_sz; ++i)
            coeff_new[i] = coeff_[i*stride_];
        if (is_owned_)
          ds::destroy_and_release(coeff_, sz_*stride_);
        coeff_ = coeff_new;
        sz_ = sz;
        stride_ = 1;
        is_owned_ = true;
      }
    }

    //! Reset storage to given array, size, and stride
    KOKKOS_INLINE_FUNCTION
    void shallowReset(pointer v, const ordinal_type& sz,
                      const ordinal_type& stride, bool owned) {
      if (is_owned_)
        ds::destroy_and_release(coeff_, sz_*stride_);
      coeff_ = v;
      sz_ = sz;
      stride_ = stride;
      is_owned_ = owned;
    }

    //! Reset storage to given array, size, and stride
    KOKKOS_INLINE_FUNCTION
    void shallowReset(pointer v, const ordinal_type& sz,
                      const ordinal_type& stride, bool owned) volatile {
      if (is_owned_)
        ds::destroy_and_release(coeff_, sz_*stride_);
      coeff_ = v;
      sz_ = sz;
      stride_ = stride;
      is_owned_ = owned;
    }

    //! Return size
    KOKKOS_INLINE_FUNCTION
    ordinal_type size() const { return sz_; }

    //! Return size
    KOKKOS_INLINE_FUNCTION
    ordinal_type size() const volatile { return sz_; }

    //! Coefficient access (avoid if possible)
    KOKKOS_INLINE_FUNCTION
    const_reference operator[] (const ordinal_type& i) const {
      return coeff_[i*stride_];
    }

    //! Coefficient access (avoid if possible)
    KOKKOS_INLINE_FUNCTION
    const_volatile_reference operator[] (const ordinal_type& i) const volatile {
      return coeff_[i*stride_];
    }

    //! Coefficient access (avoid if possible)
    KOKKOS_INLINE_FUNCTION
    reference operator[] (const ordinal_type& i) {
      return coeff_[i*stride_];
    }

    //! Coefficient access (avoid if possible)
    KOKKOS_INLINE_FUNCTION
    volatile_reference operator[] (const ordinal_type& i) volatile {
      return coeff_[i*stride_];
    }

    template <int i>
    KOKKOS_INLINE_FUNCTION
    reference getCoeff() { return coeff_[i*stride_]; }

    template <int i>
    KOKKOS_INLINE_FUNCTION
    volatile_reference getCoeff() volatile { return coeff_[i*stride_]; }

    template <int i>
    KOKKOS_INLINE_FUNCTION
    const_reference getCoeff() const { return coeff_[i*stride_]; }

    template <int i>
    KOKKOS_INLINE_FUNCTION
    const_volatile_reference getCoeff() const volatile { return coeff_[i*stride_]; }

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
    pointer coeff_;

    //! Size of array used
    ordinal_type sz_;

    //! Stride of array
    ordinal_type stride_;

    //! Do we own the array
    bool is_owned_;

  };

}

#include "Stokhos_StorageHelpers.hpp"
STOKHOS_STORAGE_HELPER_STRINGNAME_DYNAMIC(DynamicStridedStorage)

#endif // STOKHOS_DYNAMIC_STORAGE_HPP
