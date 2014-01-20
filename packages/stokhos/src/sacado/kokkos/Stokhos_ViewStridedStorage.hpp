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

#ifndef STOKHOS_VIEW_STRIDED_STORAGE_HPP
#define STOKHOS_VIEW_STRIDED_STORAGE_HPP

#include "Stokhos_DynArrayTraits.hpp"

#include "Kokkos_Macros.hpp"

#include "Sacado_Traits.hpp"
#include "Stokhos_KokkosTraits.hpp"
#include <sstream>

namespace Stokhos {

  //! Dynamic storage with view semantics and strided access
  template <typename ordinal_t, typename value_t, typename device_t>
  class ViewStridedStorage {
  public:

    static const bool is_static = false;
    static const int static_size = 0;
    static const bool supports_reset = true;

    typedef ordinal_t ordinal_type;
    typedef value_t value_type;
    typedef device_t device_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef Stokhos::DynArrayTraits<value_type,device_type> ds;

    //! Turn ViewStridedStorage into a meta-function class usable with mpl::apply
    template <typename ord_t, typename val_t = value_t , typename dev_t = device_t>
    struct apply {
      typedef ViewStridedStorage<ord_t,val_t,dev_t> type;
    };

    //! Constructor to satisfy Sacado::MP::Vector
    /*!
     * Should remove this and add disable_if's to Sacado::MP::Vector to
     * prevent accidentally creating a ViewStridedStorage object without a
     * pointer.
     */
    KOKKOS_INLINE_FUNCTION
    ViewStridedStorage(const ordinal_type& sz,
                       const value_type& x = value_type(0.0)) {}

    //! Constructor
    KOKKOS_INLINE_FUNCTION
    ViewStridedStorage(pointer v, const ordinal_type& sz, const ordinal_type& stride) :
      coeff_(v), sz_(sz), stride_(stride) {}

    //! Constructor
    KOKKOS_INLINE_FUNCTION
    ViewStridedStorage(const ViewStridedStorage& s) :
      coeff_(s.coeff_), sz_(s.sz_), stride_(s.stride_) {}

    //! Destructor
    KOKKOS_INLINE_FUNCTION
    ~ViewStridedStorage() {}

    //! Assignment operator
    KOKKOS_INLINE_FUNCTION
    ViewStridedStorage& operator=(const ViewStridedStorage& s) {
      if (&s != this) {
        coeff_ = s.coeff_;
        sz_ = s.sz_;
        stride_ = s.stride_;
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

    //! Load values to an array of values
    KOKKOS_INLINE_FUNCTION
    void load(pointer v) {
      if (stride_ == 1)
        ds::copy(coeff_, v, sz_);
      for (ordinal_type i=0; i<sz_; ++i)
        coeff_[i*stride_] = v[i];
    }

    //! Resize to new size (values are preserved)
    KOKKOS_INLINE_FUNCTION
    void resize(const ordinal_type& sz) {
      sz_ = sz;
    }

    //! Reset storage to given array, size, and stride
    KOKKOS_INLINE_FUNCTION
    void shallowReset(pointer v, const ordinal_type& sz,
                      const ordinal_type& stride, bool owned) {
      coeff_ = v;
      sz_ = sz;
      stride_ = stride;
    }

    //! Return size
    KOKKOS_INLINE_FUNCTION
    ordinal_type size() const { return sz_; }

    //! Coefficient access (avoid if possible)
    KOKKOS_INLINE_FUNCTION
    const_reference operator[] (const ordinal_type& i) const {
      return coeff_[i*stride_];
    }

    //! Coefficient access (avoid if possible)
    KOKKOS_INLINE_FUNCTION
    reference operator[] (const ordinal_type& i) {
      return coeff_[i*stride_];
    }

    template <int i>
    KOKKOS_INLINE_FUNCTION
    reference getCoeff() { return coeff_[i*stride_]; }

    template <int i>
    KOKKOS_INLINE_FUNCTION
    const_reference getCoeff() const { return coeff_[i*stride_]; }

    //! Get coefficients
    KOKKOS_INLINE_FUNCTION
    const_pointer coeff() const { return coeff_; }

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

  };

}

namespace Sacado {
  template <typename ordinal_t, typename value_t, typename device_t>
  struct StringName< Stokhos::ViewStridedStorage<ordinal_t,
                                                 value_t,
                                                 device_t> > {
    static std::string eval() {
      std::stringstream ss;
      ss << "Stokhos::ViewStridedStorage<"
         << StringName<ordinal_t>::eval() << ","
         << StringName<value_t>::eval() << ","
         << StringName<device_t>::eval() << ">";
      return ss.str();
    }
  };
}

#endif // STOKHOS_VIEW_STRIDED_STORAGE_HPP
