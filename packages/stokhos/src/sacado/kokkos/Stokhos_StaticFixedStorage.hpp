// @HEADER
// ***********************************************************************
//
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_STATIC_FIXED_STORAGE_HPP
#define STOKHOS_STATIC_FIXED_STORAGE_HPP

#include "Stokhos_StaticArrayTraits.hpp"
#include "Stokhos_MemoryTraits.hpp"

#include "Kokkos_Core_fwd.hpp"

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
    KOKKOS_INLINE_FUNCTION
    StaticFixedStorage(const ordinal_type& sz,
                       const value_type& x = value_type(0.0)) {
      ss::fill(coeff_, Num, x);
    }

    //! Constructor for creating a view (not allowed)
    KOKKOS_INLINE_FUNCTION
    StaticFixedStorage(const ordinal_type& sz, pointer v, bool owned) {}

    //! Copy constructor
    KOKKOS_INLINE_FUNCTION
    StaticFixedStorage(const StaticFixedStorage& s) {
      ss::copy(s.coeff_, coeff_, Num);
    }

    //! Copy constructor
    KOKKOS_INLINE_FUNCTION
    StaticFixedStorage(const volatile StaticFixedStorage& s) {
      ss::copy(s.coeff_, coeff_, Num);
    }

    //! Destructor
    KOKKOS_INLINE_FUNCTION
    ~StaticFixedStorage() {}

    //! Assignment operator
    KOKKOS_INLINE_FUNCTION
    StaticFixedStorage& operator=(const StaticFixedStorage& s) {
      ss::copy(s.coeff_, coeff_, Num);
      return *this;
    }

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

#if defined(KOKKOS_HAVE_CUDA)

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
    KOKKOS_INLINE_FUNCTION
    StaticFixedStorage(const ordinal_type& sz,
                       const value_type& x = value_type(0.0)) {
      ss::fill(coeff_, Num, x);
    }

    //! Constructor for creating a view (not allowed)
    KOKKOS_INLINE_FUNCTION
    StaticFixedStorage(const ordinal_type& sz, pointer v, bool owned) {}

    //! Copy constructor
    KOKKOS_INLINE_FUNCTION
    StaticFixedStorage(const StaticFixedStorage& s) {
      ss::copy(s.coeff_, coeff_, Num);
    }

    //! Copy constructor
    KOKKOS_INLINE_FUNCTION
    StaticFixedStorage(const volatile StaticFixedStorage& s) {
      ss::copy(s.coeff_, coeff_, Num);
    }

    //! Destructor
    KOKKOS_INLINE_FUNCTION
    ~StaticFixedStorage() {}

    //! Assignment operator
    KOKKOS_INLINE_FUNCTION
    StaticFixedStorage& operator=(const StaticFixedStorage& s) {
      ss::copy(s.coeff_, coeff_, Num);
      return *this;
    }

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

namespace Sacado {
  template <typename ordinal_t, typename value_t, int Num, typename device_t>
  struct StringName< Stokhos::StaticFixedStorage<ordinal_t,
                                                 value_t,
                                                 Num,
                                                 device_t> > {
    static std::string eval() {
      std::stringstream ss;
      ss << "Stokhos::StaticFixedStorage<"
         << StringName<ordinal_t>::eval() << ","
         << StringName<value_t>::eval() << ","
         << Num << ","
         << StringName<device_t>::eval() << ">";
      return ss.str();
    }
  };
}

#endif // STOKHOS_STATIC_FIXED_STORAGE_HPP
