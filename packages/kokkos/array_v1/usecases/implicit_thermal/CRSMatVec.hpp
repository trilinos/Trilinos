/** \HEADER
 *************************************************************************
 *
 *                            Kokkos
 *                 Copyright 2010 Sandia Corporation
 *
 *  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
 *  the U.S. Government retains certain rights in this software.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions are
 *  met:
 *
 *  1. Redistributions of source code must retain the above copyright
 *  notice, this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright
 *  notice, this list of conditions and the following disclaimer in the
 *  documentation and/or other materials provided with the distribution.
 *
 *  3. Neither the name of the Corporation nor the names of the
 *  contributors may be used to endorse or promote products derived from
 *  this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
 *  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
 *  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *  EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 *  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 *  PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 *  LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 *  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *************************************************************************
 */

#include <stdexcept>
#include <string>
//----------------------------------------------------------------------------

template<class Scalar, class DeviceType>
struct CRSMatVec
{
  typedef DeviceType                                       device_type;
  typedef typename device_type::size_type                  index_type ;
  typedef Scalar                                           scalar_type ;
  typedef Kokkos::MultiVectorView<index_type,  device_type> index_vector;   
  typedef Kokkos::MultiVectorView<scalar_type, device_type> scalar_vector;
  

  scalar_vector A ;
  index_vector  A_row , A_col ;
  scalar_vector x , b ;

  CRSMatVec( const scalar_vector & arg_A ,
             const index_vector  & arg_row ,
             const index_vector  & arg_col,
             const scalar_vector & arg_x ,
             const scalar_vector & arg_b )
  : A(arg_A) , A_row(arg_row) , A_col(arg_col)
  , x(arg_x) , b(arg_b)      
  {}

  KOKKOS_MACRO_DEVICE_FUNCTION
  void operator()( index_type irow ) const
  {
    const index_type jEnd = A_row(irow+1);

    double sum = 0 ;

    for ( index_type j = A_row(irow) ; j < jEnd ; ++j ) {
      sum += A(j) * x( A_col(j) );
    }

    b(irow) = sum ;
  }

  void apply() const
  {
    const index_type n = A_row.length() - 1 ;
    parallel_for( n , *this );
  }
};

//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

#if defined KOKKOS_DEVICE_CUDA

#include <cuda_runtime.h>
#include "cusparse.h"

// Specialize CRSMatVect for CUDA sparse library

template<class Scalar>
struct CRSMatVec<Scalar , Kokkos::DeviceCuda>
{
  typedef Kokkos::DeviceCuda      device_type ;
  typedef device_type::size_type  index_type ;
  typedef Kokkos::MultiVectorView<Scalar, device_type>     scalar_vector;
  typedef Kokkos::MultiVectorView<index_type, device_type> index_vector;   
  
  Scalar * value, *x , *b ;
  int * row ,* col;
  int m;
  
  cusparseStatus_t status;
  cusparseHandle_t handle;
  cusparseMatDescr_t descra;
  
  //Constructor that initalizes cusparse library and 
  //gets pointers on the device to use 
  CRSMatVec( scalar_vector & arg_value ,
             index_vector  & arg_row ,
             index_vector  & arg_col ,
             scalar_vector & arg_x ,
             scalar_vector & arg_b )
  : value( (Scalar*)arg_value.ptr_on_device() )
  , row( (int*)arg_row.ptr_on_device() )
  , col( (int*)arg_col.ptr_on_device() )
  , x( (Scalar*)arg_x.ptr_on_device() )
  , b( (Scalar*)arg_b.ptr_on_device() )
  , m( arg_row.length() - 1 )
  {
    status = cusparseCreate(&handle);
    if(status != CUSPARSE_STATUS_SUCCESS)
    {
      throw std::runtime_error( std::string("ERROR - CUSPARSE Library Initialization failed" ) );
    }

    status = cusparseCreateMatDescr(&descra);
    if(status != CUSPARSE_STATUS_SUCCESS)
    {
      throw std::runtime_error( std::string("ERROR - CUSPARSE Library Matrix descriptor failed" ) );
    }

    cusparseSetMatType(descra , CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descra , CUSPARSE_INDEX_BASE_ZERO);

  }

  void apply() const ;
};

//Compute equivalent MatVec operation using 
//cuda sparse library function

template<>
void CRSMatVec<float, Kokkos::DeviceCuda >::apply() const
{
  float alpha = 1 , beta = 0 ;
  cusparseScsrmv(handle , CUSPARSE_OPERATION_NON_TRANSPOSE , m, m , alpha , 
        descra , value , row , col , x , beta , b);
}

template<>
void CRSMatVec<double, Kokkos::DeviceCuda>::apply() const
{
  double alpha = 1 , beta = 0 ;
  cusparseDcsrmv(handle , CUSPARSE_OPERATION_NON_TRANSPOSE , m, m , alpha , 
        descra , value , row , col , x , beta , b);
}

#endif

//----------------------------------------------------------------------------

