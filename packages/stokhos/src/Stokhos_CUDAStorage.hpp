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
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
