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

#ifndef STOKHOS_HOST_LINEAR_SPARSE_3_TENSOR_HPP
#define STOKHOS_HOST_LINEAR_SPARSE_3_TENSOR_HPP

#include "KokkosArray_Host.hpp"

#include "Stokhos_Multiply.hpp"
#include "Stokhos_LinearSparse3Tensor.hpp"
#include "Stokhos_Host_TinyVec.hpp"

namespace Stokhos {

template< typename ValueType, int BlockSize >
class Multiply< LinearSparse3Tensor< ValueType , KokkosArray::Host , BlockSize > , void , void , DefaultSparseMatOps >
{
public:

  typedef KokkosArray::Host::size_type size_type ;
  typedef LinearSparse3Tensor< ValueType , KokkosArray::Host , BlockSize > tensor_type ;

  template< typename MatrixValue , typename VectorValue >
  static void apply( const tensor_type & tensor ,
                     const MatrixValue * const a ,
                     const VectorValue * const x ,
                           VectorValue * const y )
  {
    const size_type block_size = tensor_type::block_size;
    typedef TinyVec<ValueType,block_size,true> TV;
    const size_type dim = tensor.dimension();

    const ValueType c0 = tensor.value(0);
    const ValueType c1 = tensor.value(1);
    const ValueType a0 = a[0];
    const ValueType x0 = x[0];

    if (block_size > 1) {

      TV vc0(c0), vc1(c1), va0(a0), vx0(x0), vy0;
      TV ai, ai2, xi, yi;

      const MatrixValue *aa = a;
      const VectorValue *xx = x;
      VectorValue *yy = y;
      vy0.zero();

      const size_type nBlock = dim / block_size;
      const size_type iEnd = nBlock * block_size;

      if (tensor.symmetric()) {

        size_type i=0;
        for ( ; i < iEnd; i+=block_size,aa+=block_size,xx+=block_size,yy+=block_size) {
          ai.aligned_load(aa);
          ai2 = ai;
          xi.aligned_load(xx);
          yi.aligned_load(yy);

          // y[i] += c1*(a0*xi + ai*x0);
          ai.times_equal(vx0);
          ai2.times_equal(xi);
          xi.times_equal(va0);
          xi.plus_equal(ai);
          xi.times_equal(vc1);
          yi.plus_equal(xi);
          yi.aligned_scatter(yy);

          // y0  += c1*ai*xi;
          ai2.times_equal(vc1);
          vy0.plus_equal(ai2);
        }
        ValueType y0 = vy0.sum();

        // Do remaining entries with a scalar loop
        for ( ; i < dim; ++i) {
          const ValueType ai = *aa++;
          const ValueType xi = *xx++;
          *yy++ += c1*(a0*xi + ai*x0);
          y0  += c1*ai*xi;
        }
        y[0] += y0 - 2*a0*x0;
      }
      else {

        const ValueType c2 = tensor.value(2);
        TV vc2(c2);
        size_type i=0;
        for ( ; i < iEnd; i+=block_size,aa+=block_size,xx+=block_size,yy+=block_size) {
          ai.aligned_load(aa);
          ai2 = ai;
          xi.aligned_load(xx);
          yi.aligned_load(yy);

          // y[i] += c1*(a0*xi + ai*x0) + c2*aixi;
          ai.times_equal(vx0);
          ai2.times_equal(xi);
          xi.times_equal(va0);
          xi.plus_equal(ai);
          xi.times_equal(vc1);
          yi.plus_equal(xi);
          ai = ai2;
          ai.times_equal(vc2);
          yi.plus_equal(ai);
          yi.aligned_scatter(yy);

          // y0  += c1*aixi;
          ai2.times_equal(vc2);
          vy0.plus_equal(ai2);
        }
        ValueType y0 = vy0.sum();

        // Do remaining entries with a scalar loop
        for ( ; i < dim; ++i) {
          const ValueType ai = *aa++;
          const ValueType xi = *xx++;
          const ValueType aixi = ai*xi;
          *yy++ += c1*(a0*xi + ai*x0) + c2*aixi;
          y0  += c1*aixi;
        }
        y[0] += y0 - 2*a0*x0;

      }

    }

    else {

      ValueType y0 = c0*a0*x0;

      if (tensor.symmetric()) {

        for ( size_type i = 1; i < dim; ++i) {
          const ValueType ai = a[i];
          const ValueType xi = x[i];
          y[i] += c1*(a0*xi + ai*x0);
          y0  += c1*ai*xi;
        }
        y[0] += y0;

      }
      else {

        const ValueType c2 = tensor.value(2);
        for ( size_type i = 1; i < dim; ++i) {
          const ValueType ai = a[i];
          const ValueType xi = x[i];
          const ValueType aixi = ai*xi;
          y[i] += c1*(a0*xi + ai*x0) + c2*aixi;
          y0  += c1*aixi;
        }
        y[0] += y0;

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

#endif /* #ifndef STOKHOS_HOST_LINEAR_SPARSE_3_TENSOR_HPP */
