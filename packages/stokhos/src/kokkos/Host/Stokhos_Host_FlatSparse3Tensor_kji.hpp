// @HEADER
// ***********************************************************************
// 
//                     Stokhos Package
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
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef STOKHOS_HOST_FLAT_SPARSE_3_TENSOR_KJI_HPP
#define STOKHOS_HOST_FLAT_SPARSE_3_TENSOR_KJI_HPP

#include "KokkosArray_Host.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_FlatSparse3Tensor_kji.hpp"

namespace Stokhos {

template< typename ValueType >
class Multiply< FlatSparse3Tensor_kji< ValueType , KokkosArray::Host > , void , void >
{
public:
  
  typedef KokkosArray::Host::size_type size_type ;
  typedef FlatSparse3Tensor_kji< ValueType , KokkosArray::Host > tensor_type ;

  template< typename MatrixValue , typename VectorValue >
  static void apply( const tensor_type & tensor ,
                     const MatrixValue * const a ,
                     const VectorValue * const x ,
                           VectorValue * const y )
  {
    const size_type nk = tensor.num_k();

    // Loop over k
    for ( size_type k = 0; k < nk; ++k) {
      const MatrixValue ak = a[k];
      const VectorValue xk = x[k];

      // Loop over j for this k
      const size_type nj = tensor.num_j(k);
      const size_type jBeg = tensor.j_begin(k);
      const size_type jEnd = jBeg + nj;
      for (size_type jEntry = jBeg; jEntry < jEnd; ++jEntry) {
	const size_type j = tensor.j_coord(jEntry);
	VectorValue tmp = a[j] * xk + ak * x[j];

	// Loop over i for this k,j
	const size_type ni = tensor.num_i(jEntry);
	const size_type iBeg = tensor.i_begin(jEntry);
	const size_type iEnd = iBeg + ni;
	for (size_type iEntry = iBeg; iEntry < iEnd; ++iEntry) {
	  const size_type i = tensor.i_coord(iEntry);
	  y[i] += tensor.value(iEntry) * tmp;
	}
      }
    }
  }

  static size_type matrix_size( const tensor_type & tensor )
  { return tensor.dimension(); }

  static size_type vector_size( const tensor_type & tensor )
  { return tensor.dimension(); }
};

//----------------------------------------------------------------------------

} // namespace Stokhos

#endif /* #ifndef STOKHOS_HOST_SPARSEPRODUCTTENSOR_KJI_HPP */

