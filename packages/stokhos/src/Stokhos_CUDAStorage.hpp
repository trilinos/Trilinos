// $Id$
// $Source$ 
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

#ifndef STOKHOS_CUDA_STORAGE_HPP
#define STOKHOS_CUDA_STORAGE_HPP

#include "Stokhos_ConfigDefs.h"

#ifdef HAVE_STOKHOS_THRUST

#include <thrust/device_vector.h>

namespace Stokhos {

  template <typename ordinal_type, typename value_type>
  class CUDAStorage {
  public:

    typedef typename thrust::device_vector<value_type>::reference reference;
    typedef typename thrust::device_vector<value_type>::const_reference const_reference;
    typedef typename thrust::device_vector<value_type>::pointer pointer;
    typedef typename thrust::device_vector<value_type>::const_pointer const_pointer;

    //! Turn CUDAStorage into a meta-function class usable with mpl::apply
    template <typename ord_t, typename val_t> 
    struct apply {
      typedef CUDAStorage<ord_t, val_t> type;
    };
    
    //! Constructor
    CUDAStorage(const ordinal_type& sz);

    //! Copy constructor
    CUDAStorage(const CUDAStorage& s);

    //! Assignment operator
    CUDAStorage& operator=(const CUDAStorage& s);

    //! Destructor
    ~CUDAStorage();

    //! Initialize values to a constant value
    void init(const value_type& v);

    //! Initialize values to an array of values
    void init(const value_type* v, const ordinal_type& sz = 0);

    //! Load coefficients to an array of values
    void load(value_type* v);

    //! Resize to new size (values are preserved)
    void resize(const ordinal_type& sz);

    //! Return size
    ordinal_type size() const;

    //! Coefficient access (avoid if possible)
    const_reference operator[] (const ordinal_type& i) const;

    //! Coefficient access (avoid if possible)
    reference operator[] (const ordinal_type& i);

    //! Get coefficients
    const_pointer coeff() const;

    //! Get coefficients
    pointer coeff();

    //! Get coefficient device vector
    thrust::device_vector<value_type>& dev_vector();

    //! Get coefficient device vector
    const thrust::device_vector<value_type>& dev_vector() const;

  private:

    //! Coefficient values
    thrust::device_vector<value_type> coeff_;

  };

}

#endif // HAVE_STOKHOS_THRUST

#endif // STOKHOS_CUDA_STORAGE_HPP
