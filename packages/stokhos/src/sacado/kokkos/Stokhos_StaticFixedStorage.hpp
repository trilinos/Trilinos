// @HEADER
// *****************************************************************************
//                           Stokhos Package
//
// Copyright 2009 NTESS and the Stokhos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_STATIC_FIXED_STORAGE_HPP
#define STOKHOS_STATIC_FIXED_STORAGE_HPP

#include "Stokhos_StaticArrayTraits.hpp"
#include "Stokhos_MemoryTraits.hpp"

// We are hooking into Kokkos Core internals here
// Need to define this macro since we include non-public headers
#ifndef KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE
#define KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#endif
#include "Kokkos_Core_fwd.hpp"
#include "Cuda/Kokkos_Cuda.hpp"
#ifdef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE
#undef KOKKOS_IMPL_PUBLIC_INCLUDE_NOTDEFINED_CORE
#endif

#include "Sacado_Traits.hpp"
#include "Stokhos_KokkosTraits.hpp"
#include <sstream>

namespace Stokhos {

  //! Statically allocated storage class
  template <typename ordinal_t, typename value_t, int Num, typename device_t>
  class StaticFixedStorage {
  public:

    static const bool is_static = true ;
    static const int static_size = Num;
    static const bool supports_reset = false;

    typedef ordinal_t ordinal_type;
    typedef value_t value_type;
    typedef typename device_t::execution_space execution_space;
    typedef typename device_t::memory_space memory_space;
    typedef value_type& reference;
    typedef volatile value_type& volatile_reference;
    typedef const value_type& const_reference;
    typedef const volatile value_type& const_volatile_reference;
    typedef value_type* pointer;
    typedef volatile value_type* volatile_pointer;
    typedef const value_type* const_pointer;
    typedef const volatile value_type* const_volatile_pointer;
    typedef Stokhos::StaticArrayTraits<value_type,execution_space> ss;

    typedef typename Stokhos::MemoryTraits<memory_space> MemTraits;

    //! Turn StaticFixedStorage into a meta-function class usable with mpl::apply
    template <typename ord_t, typename val_t = value_t , typename dev_t = device_t >
    struct apply {
      typedef StaticFixedStorage<ord_t,val_t,Num,dev_t> type;
    };

    template <int N>
    struct apply_N {
      typedef StaticFixedStorage<ordinal_type,value_type,N,device_t> type;
    };

    //! Constructor
    KOKKOS_DEFAULTED_FUNCTION
    StaticFixedStorage() = default;

    //! Constructor
    KOKKOS_INLINE_FUNCTION
    StaticFixedStorage(const ordinal_type& sz,
                       const value_type& x = value_type(0.0)) {
      ss::fill(coeff_, Num, x);
    }

    //! Constructor from array
    KOKKOS_INLINE_FUNCTION
    StaticFixedStorage(const ordinal_type& sz, const value_type* x) {
      ss::copy(x, coeff_, sz);
    }

    //! Constructor for creating a view (not allowed)
    KOKKOS_INLINE_FUNCTION
    StaticFixedStorage(const ordinal_type& sz, pointer v, bool owned) {}

    //! Copy constructor
    KOKKOS_DEFAULTED_FUNCTION
    StaticFixedStorage(const StaticFixedStorage& s) = default;

    //! Copy constructor
    KOKKOS_INLINE_FUNCTION
    StaticFixedStorage(const volatile StaticFixedStorage& s) {
      ss::copy(s.coeff_, coeff_, Num);
    }

    //! Destructor
    KOKKOS_DEFAULTED_FUNCTION
    ~StaticFixedStorage() = default;

    //! Assignment operator
    KOKKOS_DEFAULTED_FUNCTION
    StaticFixedStorage& operator=(const StaticFixedStorage& s) = default;

    //! Assignment operator
    KOKKOS_INLINE_FUNCTION
    StaticFixedStorage& operator=(const volatile StaticFixedStorage& s) {
      ss::copy(s.coeff_, coeff_, Num);
      return *this;
    }

    //! Assignment operator
    KOKKOS_INLINE_FUNCTION
    /*volatile*/ StaticFixedStorage&
    operator=(const StaticFixedStorage& s) volatile {
      ss::copy(s.coeff_, coeff_, Num);
      return const_cast<StaticFixedStorage&>(*this);
    }

    //! Assignment operator
    KOKKOS_INLINE_FUNCTION
    /*volatile*/ StaticFixedStorage&
    operator=(const volatile StaticFixedStorage& s) volatile {
      ss::copy(s.coeff_, coeff_, Num);
      return const_cast<StaticFixedStorage&>(*this);
    }

    //! Initialize values to a constant value
    KOKKOS_INLINE_FUNCTION
    void init(const_reference v) {
      ss::fill(coeff_, Num, v);
    }

    //! Initialize values to a constant value
    KOKKOS_INLINE_FUNCTION
    void init(const_reference v) volatile {
      ss::fill(coeff_, Num, v);
    }

    //! Initialize values to an array of values
    KOKKOS_INLINE_FUNCTION
    void init(const_pointer v, const ordinal_type& sz = 0) {
      if (sz == 0)
        ss::copy(v, coeff_, Num);
      else
        ss::copy(v, coeff_, sz);
    }

    //! Initialize values to an array of values
    KOKKOS_INLINE_FUNCTION
    void init(const_pointer v, const ordinal_type& sz = 0) volatile {
      if (sz == 0)
        ss::copy(v, coeff_, Num);
      else
        ss::copy(v, coeff_, sz);
    }

    //! Load values to an array of values
    KOKKOS_INLINE_FUNCTION
    void load(pointer v) {
      ss::copy(v, coeff_, Num);
    }

    //! Load values to an array of values
    KOKKOS_INLINE_FUNCTION
    void load(pointer v) volatile {
      ss::copy(v, coeff_, Num);
    }

    //! Resize to new size (values are preserved)
    KOKKOS_INLINE_FUNCTION
    void resize(const ordinal_type& sz) {}

    //! Resize to new size (values are preserved)
    KOKKOS_INLINE_FUNCTION
    void resize(const ordinal_type& sz) volatile {}

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
    static ordinal_type size() { return Num; }

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
      return coeff_[i]; }

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
#if STOKHOS_ALIGN_MEMORY && ( defined(__INTEL_COMPILER) )
    value_type coeff_[Num] __attribute__((aligned(MemTraits::Alignment)));
#else
    value_type coeff_[Num];
#endif

  };

#if defined(KOKKOS_ENABLE_CUDA)

  //! Statically allocated storage class
  /*!
   * Partial specialization for Cuda where there is no alignment specified
   * for the coeff_ array.
   */
  template <typename ordinal_t, typename value_t, int Num>
  class StaticFixedStorage<ordinal_t,value_t,Num,Kokkos::Cuda> {
  public:

    static const bool is_static = true ;
    static const int static_size = Num;
    static const bool supports_reset = false;

    typedef ordinal_t ordinal_type;
    typedef value_t value_type;
    typedef Kokkos::Cuda execution_space;
    typedef value_type& reference;
    typedef volatile value_type& volatile_reference;
    typedef const value_type& const_reference;
    typedef const volatile value_type& const_volatile_reference;
    typedef value_type* pointer;
    typedef volatile value_type* volatile_pointer;
    typedef const value_type* const_pointer;
    typedef const volatile value_type* const_volatile_pointer;
    typedef Stokhos::StaticArrayTraits<value_type,execution_space> ss;

    typedef typename execution_space::memory_space memory_space;
    typedef typename Stokhos::MemoryTraits<memory_space> MemTraits;

    //! Turn StaticFixedStorage into a meta-function class usable with mpl::apply
    template <typename ord_t, typename val_t = value_t , typename dev_t = execution_space >
    struct apply {
      typedef StaticFixedStorage<ord_t,val_t,Num,dev_t> type;
    };

    template <int N>
    struct apply_N {
      typedef StaticFixedStorage<ordinal_type,value_type,N,execution_space> type;
    };

    //! Constructor
    KOKKOS_DEFAULTED_FUNCTION
    StaticFixedStorage() = default;

    //! Constructor
    KOKKOS_INLINE_FUNCTION
    StaticFixedStorage(const ordinal_type& sz,
                       const value_type& x = value_type(0.0)) {
      ss::fill(coeff_, Num, x);
    }

    //! Constructor from array
    KOKKOS_INLINE_FUNCTION
    StaticFixedStorage(const ordinal_type& sz, const value_type* x) {
      ss::copy(x, coeff_, sz);
    }

    //! Constructor for creating a view (not allowed)
    KOKKOS_INLINE_FUNCTION
    StaticFixedStorage(const ordinal_type& sz, pointer v, bool owned) {}

    //! Copy constructor
    KOKKOS_DEFAULTED_FUNCTION
    StaticFixedStorage(const StaticFixedStorage& s) = default;

    //! Copy constructor
    KOKKOS_INLINE_FUNCTION
    StaticFixedStorage(const volatile StaticFixedStorage& s) {
      ss::copy(s.coeff_, coeff_, Num);
    }

    //! Destructor
    KOKKOS_DEFAULTED_FUNCTION
    ~StaticFixedStorage() = default;

    //! Assignment operator
    KOKKOS_DEFAULTED_FUNCTION
    StaticFixedStorage& operator=(const StaticFixedStorage& s) = default;

    //! Assignment operator
    KOKKOS_INLINE_FUNCTION
    StaticFixedStorage& operator=(const volatile StaticFixedStorage& s) {
      ss::copy(s.coeff_, coeff_, Num);
      return *this;
    }

    //! Assignment operator
    KOKKOS_INLINE_FUNCTION
    /*volatile*/ StaticFixedStorage&
    operator=(const StaticFixedStorage& s) volatile {
      ss::copy(s.coeff_, coeff_, Num);
      return const_cast<StaticFixedStorage&>(*this);
    }

    //! Assignment operator
    KOKKOS_INLINE_FUNCTION
    /*volatile*/ StaticFixedStorage&
    operator=(const volatile StaticFixedStorage& s) volatile {
      ss::copy(s.coeff_, coeff_, Num);
      return const_cast<StaticFixedStorage&>(*this);
    }

    //! Initialize values to a constant value
    KOKKOS_INLINE_FUNCTION
    void init(const_reference v) {
      ss::fill(coeff_, Num, v);
    }

    //! Initialize values to a constant value
    KOKKOS_INLINE_FUNCTION
    void init(const_reference v) volatile {
      ss::fill(coeff_, Num, v);
    }

    //! Initialize values to an array of values
    KOKKOS_INLINE_FUNCTION
    void init(const_pointer v, const ordinal_type& sz = 0) {
      if (sz == 0)
        ss::copy(v, coeff_, Num);
      else
        ss::copy(v, coeff_, sz);
    }

    //! Initialize values to an array of values
    KOKKOS_INLINE_FUNCTION
    void init(const_pointer v, const ordinal_type& sz = 0) volatile {
      if (sz == 0)
        ss::copy(v, coeff_, Num);
      else
        ss::copy(v, coeff_, sz);
    }

    //! Load values to an array of values
    KOKKOS_INLINE_FUNCTION
    void load(pointer v) {
      ss::copy(v, coeff_, Num);
    }

    //! Load values to an array of values
    KOKKOS_INLINE_FUNCTION
    void load(pointer v) volatile {
      ss::copy(v, coeff_, Num);
    }

    //! Resize to new size (values are preserved)
    KOKKOS_INLINE_FUNCTION
    void resize(const ordinal_type& sz) {}

    //! Resize to new size (values are preserved)
    KOKKOS_INLINE_FUNCTION
    void resize(const ordinal_type& sz) volatile {}

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
    static ordinal_type size() { return Num; }

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
      return coeff_[i]; }

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

  };

#endif

}

#include "Stokhos_StorageHelpers.hpp"
STOKHOS_STORAGE_HELPER_STRINGNAME_STATIC(StaticFixedStorage)

#endif // STOKHOS_STATIC_FIXED_STORAGE_HPP
