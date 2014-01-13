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

#ifndef STOKHOS_DYNAMIC_THREADED_STORAGE_HPP
#define STOKHOS_DYNAMIC_THREADED_STORAGE_HPP

#include "Stokhos_DynArrayTraits.hpp"

#include "Sacado_Traits.hpp"
#include "Stokhos_KokkosTraits.hpp"
#include <sstream>

namespace Stokhos {

  //! Dynamically allocated storage class with striding
  template <typename ordinal_t, typename value_t, typename device_t>
  class DynamicThreadedStorage {
  public:

    typedef ordinal_t ordinal_type;
    typedef value_t value_type;
    typedef device_t device_type;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef Stokhos::DynArrayTraits<value_type,device_type> ds;

    //! Turn DynamicThreadedStorage into a meta-function class usable with mpl::apply
    template <typename ord_t, typename val_t = value_t, typename dev_t = device_t>
    struct apply {
      typedef DynamicThreadedStorage<ord_t,val_t,dev_t> type;
    };

    //! Constructor
    DynamicThreadedStorage(const ordinal_type& sz,
                  const value_type& x = value_type(0.0));

    //! Copy constructor
    DynamicThreadedStorage(const DynamicThreadedStorage& s);

    //! Destructor
    ~DynamicThreadedStorage();

    //! Assignment operator
    DynamicThreadedStorage& operator=(const DynamicThreadedStorage& s);

    //! Initialize values to a constant value
    void init(const_reference v);

    //! Initialize values to an array of values
    void init(const_pointer v, const ordinal_type& sz_ = 0);

    //! Load values to an array of values
    void load(pointer v);

    //! Resize to new size (values are preserved)
    void resize(const ordinal_type& sz);

    //! Reset storage to given array, size, and stride
    void shallowReset(pointer v, const ordinal_type& sz,
                      const ordinal_type& stride, bool owned);

    //! Return size
    static ordinal_type size();

    //! Coefficient access (avoid if possible)
    const_reference operator[] (const ordinal_type& i) const;

    //! Coefficient access (avoid if possible)
    reference operator[] (const ordinal_type& i);

    //! Get coefficients
    const_pointer coeff() const;

    //! Get coefficients
    pointer coeff();

  };

}

namespace Sacado {
  template <typename ordinal_t, typename value_t, typename device_t>
  struct StringName< Stokhos::DynamicThreadedStorage<ordinal_t,
                                                     value_t,
                                                     device_t> > {
    static std::string eval() {
      std::stringstream ss;
      ss << "Stokhos::DynamicThreadedStorage<"
         << StringName<ordinal_t>::eval() << ","
         << StringName<value_t>::eval() << ","
         << StringName<device_t>::eval() << ">";
      return ss.str();
    }
  };
}

// No Host specialization

// Cuda specialization
#include "Kokkos_Cuda.hpp"
#include "Stokhos_DynamicThreadedStorage_cuda.hpp"

#endif // STOKHOS_DYNAMIC_STORAGE_HPP
