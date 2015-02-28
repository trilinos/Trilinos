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

#ifndef STOKHOS_DYNAMIC_STORAGE_HPP
#define STOKHOS_DYNAMIC_STORAGE_HPP

#include "Stokhos_DynArrayTraits.hpp"

#include "Kokkos_Macros.hpp"

#include "Sacado_Traits.hpp"
#include "Stokhos_KokkosTraits.hpp"
#include <sstream>

namespace Stokhos {

  template <typename ordinal_t, typename value_t, typename device_t>
  class DynamicStorage {
  public:

    static const bool is_static = false;
    static const int static_size = 0;
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
    typedef Stokhos::DynArrayTraits<value_type,execution_space> ds;

    //! Turn DynamicStorage into a meta-function class usable with mpl::apply
    template <typename ord_t, typename val_t = value_t , typename dev_t = device_t >
    struct apply {
      typedef DynamicStorage<ord_t,val_t,dev_t> type;
    };

    template <int N>
    struct apply_N {
      typedef DynamicStorage<ordinal_type,value_type,device_t> type;
    };

    //! Constructor
    KOKKOS_INLINE_FUNCTION
    DynamicStorage(const ordinal_type& sz,
                   const value_type& x = value_type(0.0)) :
      sz_(sz), is_view_(false) {
      if (sz_ > 1) {
        coeff_ = ds::get_and_fill(sz_, x);
        is_constant_ = false;
      }
      else {
        coeff_0_ = x;
        coeff_ = &coeff_0_;
        is_constant_ = true;
      }
    }

    //! Constructor for creating a view
    // This might be problematic if sz == 1 and owned == true
    KOKKOS_INLINE_FUNCTION
    DynamicStorage(const ordinal_type& sz, pointer v, bool owned) :
      coeff_(v), sz_(sz), is_view_(!owned), is_constant_(false) {}

    //! Constructor
    KOKKOS_INLINE_FUNCTION
    DynamicStorage(const DynamicStorage& s) :
      sz_(s.sz_), is_view_(false) {
      if (sz_ > 1) {
        coeff_ = ds::get_and_fill(s.coeff_, sz_);
        is_constant_ = false;
      }
      else {
        coeff_0_ = s.coeff_[0];
        coeff_ = &coeff_0_;
        is_constant_ = true;
      }
    }

    //! Constructor
    KOKKOS_INLINE_FUNCTION
    DynamicStorage(const volatile DynamicStorage& s) :
      sz_(s.sz_), is_view_(false) {
      if (sz_ > 1) {
        coeff_ = ds::get_and_fill(s.coeff_, sz_);
        is_constant_ = false;
      }
      else {
        coeff_0_ = s.coeff_[0];
        coeff_ = &coeff_0_;
        is_constant_ = true;
      }
    }

    //! Destructor
    KOKKOS_INLINE_FUNCTION
    ~DynamicStorage() {
      if (!is_view_ && !is_constant_) ds::destroy_and_release(coeff_, sz_);
    }

    //! Assignment operator
    // To do:  add error check if is_view_ == true && s.sz_ > sz_
    KOKKOS_INLINE_FUNCTION
    DynamicStorage& operator=(const DynamicStorage& s) {
      if (&s != this) {
        // Only reallocate if we own the array and the sizes
        // differ
        if (!is_view_ && s.sz_ != sz_) {
          if (!is_constant_)
            ds::destroy_and_release(coeff_, sz_);
          if (s.sz_ > 1) {
            coeff_ = ds::get_and_fill(s.coeff_, s.sz_);
            is_constant_ = false;
          }
          else {
            coeff_0_ = s.coeff_[0];
            coeff_ = &coeff_0_;
            is_constant_ = true;
          }
          sz_ = s.sz_;
        }
        else {
          // This should work if is_view_ or s.sz_ == sz_, regardless
          // of whether s.sz_ or sz_ is 1
          ds::copy(s.coeff_, coeff_, s.sz_);
        }
      }
      return *this;
    }

    //! Assignment operator
    // To do:  add error check if is_view_ == true && s.sz_ > sz_
    KOKKOS_INLINE_FUNCTION
    DynamicStorage& operator=(const volatile DynamicStorage& s) {
      if (&s != this) {
        // Only reallocate if we own the array and the sizes
        // differ
        if (!is_view_ && s.sz_ != sz_) {
          if (!is_constant_)
            ds::destroy_and_release(coeff_, sz_);
          if (s.sz_ > 1) {
            coeff_ = ds::get_and_fill(s.coeff_, s.sz_);
            is_constant_ = false;
          }
          else {
            coeff_0_ = s.coeff_[0];
            coeff_ = &coeff_0_;
            is_constant_ = true;
          }
          sz_ = s.sz_;
        }
        else {
          // This should work if is_view_ or s.sz_ == sz_, regardless
          // of whether s.sz_ or sz_ is 1
          ds::copy(s.coeff_, coeff_, s.sz_);
        }
      }
      return *this;
    }

    //! Assignment operator
    // To do:  add error check if is_view_ == true && s.sz_ > sz_
    KOKKOS_INLINE_FUNCTION
    /*volatile*/ DynamicStorage& operator=(const DynamicStorage& s) volatile {
      if (&s != this) {
        // Only reallocate if we own the array and the sizes
        // differ
        if (!is_view_ && s.sz_ != sz_) {
          if (!is_constant_)
            ds::destroy_and_release(coeff_, sz_);
          if (s.sz_ > 1) {
            coeff_ = ds::get_and_fill(s.coeff_, s.sz_);
            is_constant_ = false;
          }
          else {
            coeff_0_ = s.coeff_[0];
            coeff_ = &coeff_0_;
            is_constant_ = true;
          }
          sz_ = s.sz_;
        }
        else {
          // This should work if is_view_ or s.sz_ == sz_, regardless
          // of whether s.sz_ or sz_ is 1
          ds::copy(s.coeff_, coeff_, s.sz_);
        }
      }
      return const_cast<DynamicStorage&>(*this);
    }

    //! Assignment operator
    // To do:  add error check if is_view_ == true && s.sz_ > sz_
    KOKKOS_INLINE_FUNCTION
    /*volatile*/ DynamicStorage&
    operator=(const volatile DynamicStorage& s) volatile {
      if (&s != this) {
        // Only reallocate if we own the array and the sizes
        // differ
        if (!is_view_ && s.sz_ != sz_) {
          if (!is_constant_)
            ds::destroy_and_release(coeff_, sz_);
          if (s.sz_ > 1) {
            coeff_ = ds::get_and_fill(s.coeff_, s.sz_);
            is_constant_ = false;
          }
          else {
            coeff_0_ = s.coeff_[0];
            coeff_ = &coeff_0_;
            is_constant_ = true;
          }
          sz_ = s.sz_;
        }
        else {
          // This should work if is_view_ or s.sz_ == sz_, regardless
          // of whether s.sz_ or sz_ is 1
          ds::copy(s.coeff_, coeff_, s.sz_);
        }
      }
      return const_cast<DynamicStorage&>(*this);
    }

    //! Initialize values to a constant value
    KOKKOS_INLINE_FUNCTION
    void init(const_reference v) {
      ds::fill(coeff_, sz_, v);
    }

    //! Initialize values to a constant value
    KOKKOS_INLINE_FUNCTION
    void init(const_reference v) volatile {
      ds::fill(coeff_, sz_, v);
    }

    //! Initialize values to an array of values
    KOKKOS_INLINE_FUNCTION
    void init(const_pointer v, const ordinal_type& sz = 0) {
      if (sz == 0)
        ds::copy(v, coeff_, sz_);
      else
        ds::copy(v, coeff_, sz);
    }

    //! Initialize values to an array of values
    KOKKOS_INLINE_FUNCTION
    void init(const_pointer v, const ordinal_type& sz = 0) volatile {
      if (sz == 0)
        ds::copy(v, coeff_, sz_);
      else
        ds::copy(v, coeff_, sz);
    }

    //! Load values to an array of values
    KOKKOS_INLINE_FUNCTION
    void load(pointer v) {
      ds::copy(v, coeff_, sz_);
    }

    //! Load values to an array of values
    KOKKOS_INLINE_FUNCTION
    void load(pointer v) volatile {
      ds::copy(v, coeff_, sz_);
    }

    //! Resize to new size (values are preserved)
    KOKKOS_INLINE_FUNCTION
    void resize(const ordinal_type& sz) {
      if (!is_view_ && sz != sz_) {
        if (sz > 1) {
          value_type *coeff_new = ds::get_and_fill(sz);
          if (sz > sz_)
            ds::copy(coeff_, coeff_new, sz_);
          else
            ds::copy(coeff_, coeff_new, sz);
          if (!is_view_ && !is_constant_)
            ds::destroy_and_release(coeff_, sz_);
          coeff_ = coeff_new;
          is_constant_ = false;
        }
        else {
          coeff_0_ = coeff_[0];
          if (!is_constant_)
            ds::destroy_and_release(coeff_, sz_);
          coeff_ = &coeff_0_;
          is_constant_ = true;
        }
        sz_ = sz;
      }
    }

    //! Resize to new size (values are preserved)
    KOKKOS_INLINE_FUNCTION
    void resize(const ordinal_type& sz) volatile {
      if (!is_view_ && sz != sz_) {
        if (sz > 1) {
          value_type *coeff_new = ds::get_and_fill(sz);
          if (sz > sz_)
            ds::copy(coeff_, coeff_new, sz_);
          else
            ds::copy(coeff_, coeff_new, sz);
          if (!is_view_ && !is_constant_)
            ds::destroy_and_release(coeff_, sz_);
          coeff_ = coeff_new;
          is_constant_ = false;
        }
        else {
          coeff_0_ = coeff_[0];
          if (!is_constant_)
            ds::destroy_and_release(coeff_, sz_);
          coeff_ = &coeff_0_;
          is_constant_ = true;
        }
        sz_ = sz;
      }
    }

    //! Reset storage to given array, size, and stride
    KOKKOS_INLINE_FUNCTION
    void shallowReset(pointer v, const ordinal_type& sz,
                      const ordinal_type& stride, bool owned) {
      if (!is_view_ && !is_constant_)
        ds::destroy_and_release(coeff_, sz_);
      coeff_ = v;
      sz_ = sz;
      is_view_ = !owned;
    }

    //! Reset storage to given array, size, and stride
    KOKKOS_INLINE_FUNCTION
    void shallowReset(pointer v, const ordinal_type& sz,
                      const ordinal_type& stride, bool owned) volatile {
      if (!is_view_ && !is_constant_)
        ds::destroy_and_release(coeff_, sz_);
      coeff_ = v;
      sz_ = sz;
      is_view_ = !owned;
    }

    //! Return size
    KOKKOS_INLINE_FUNCTION
    ordinal_type size() const { return sz_; }

    //! Return size
    KOKKOS_INLINE_FUNCTION
    ordinal_type size() const volatile { return sz_; }

    //! Return whether storage is a view
    KOKKOS_INLINE_FUNCTION
    bool is_view() const { return is_view_; }

     //! Return whether storage is a view
    KOKKOS_INLINE_FUNCTION
    bool is_view() const volatile { return is_view_; }

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
    pointer coeff_;

    //! Coefficient value when sz_ == 1 (i.e., a constant)
    /*! Eliminates dynamic memory allocation for this common use-case */
    value_type coeff_0_;

    //! Size of array used
    ordinal_type sz_;

    //! Do we own the array
    bool is_view_;

    //! Is the coefficient array length-1
    bool is_constant_;

  };

}

namespace Sacado {
  template <typename ordinal_t, typename value_t, typename device_t>
  struct StringName< Stokhos::DynamicStorage<ordinal_t,
                                             value_t,
                                             device_t> > {
    static std::string eval() {
      std::stringstream ss;
      ss << "Stokhos::DynamicStorage<"
         << StringName<ordinal_t>::eval() << ","
         << StringName<value_t>::eval() << ","
         << StringName<device_t>::eval() << ">";
      return ss.str();
    }
  };
}

#endif // STOKHOS_DYNAMIC_STORAGE_HPP
