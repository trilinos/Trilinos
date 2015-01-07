/*
//@HEADER
// ************************************************************************
//
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#ifndef KOKKOS_CRSMATRIX_MKL_HPP
#define KOKKOS_CRSMATRIX_MKL_HPP

/// \file Kokkos_CrsMatrix_MKL.hpp
/// \brief Intel MKL version of Kokkos' sparse matrix interface

#include <impl/Kokkos_PhysicalLayout.hpp>

namespace Kokkos {

/// \brief Attempt to compute y = s_a * A * x using Intel MKL.
///
/// If the given combination of scalar types is not supported, do
/// nothing and return false; otherwise, compute y = s_a * A * x and
/// return true.
template<typename T, class RangeVector,class CrsMatrix,class DomainVector>
bool
MV_Multiply_DoMKL (typename Kokkos::Impl::enable_if<
		     ! Kokkos::Impl::is_same<T,double>::value && 
		     ! Kokkos::Impl::is_same<T,float>::value, 
		     typename RangeVector::value_type>::type s_b,
		   const RangeVector& y, 
		   typename DomainVector::value_type s_a,
		   const CrsMatrix& A, 
		   const DomainVector& x) 
{
  return false;
}

template<typename T, class RangeVector,class CrsMatrix,class DomainVector>
bool MV_Multiply_DoMKL(typename Kokkos::Impl::enable_if<Kokkos::Impl::is_same<T,double>::value, double  >::type s_b
                ,const RangeVector & y, double s_a,
                const CrsMatrix & A , const DomainVector & x) {

  char matdescra[6] = "GLNC0";
  char transa = 'N';
  int m = A.numRows();
  int n = x.dimension_1();
  int k = A.numCols();
  double* x_ptr = (double*)x.ptr_on_device();
  double* y_ptr = (double*)y.ptr_on_device();
  if(x.dimension_1()>1) {
    Impl::PhysicalLayout layout_x(x);
    Impl::PhysicalLayout layout_y(y);
    if((layout_x.layout_type!=layout_x.Right) || layout_y.layout_type!=layout_y.Right) return false;

    int stride_x = layout_x.stride[0];
    int stride_y = layout_y.stride[0];
    // FIXME (mfh 09 Aug 2013) Doesn't this interface only work with
    // row-major multivectors?  I recall that only the "Fortran"
    // version works with column-major multivectors, and it requires
    // that the column indices be one-based.  See
    // KokkosClassic::MklSparseOps for an example.
    mkl_dcsrmm(&transa,
               &m, &n, &k,
               &s_a,
               matdescra,
               A.values.ptr_on_device(),
               A.graph.entries.ptr_on_device(),
               (int*) &A.graph.row_map(0),
               (int*) &A.graph.row_map(1),
               x_ptr,
               &stride_x,
               &s_b,
               y_ptr,
               &stride_y);
    } else
      mkl_dcsrmv(&transa,
                 &m, &k,
                 &s_a,
                 matdescra,
                 A.values.ptr_on_device(),
                 A.graph.entries.ptr_on_device(),
                 (int*) &A.graph.row_map(0),
                 (int*) &A.graph.row_map(1),
                 x_ptr,
                 &s_b,
                 y_ptr);
  return true;
}

template<typename T, class RangeVector,class CrsMatrix,class DomainVector>
bool MV_Multiply_DoMKL(typename Kokkos::Impl::enable_if<Kokkos::Impl::is_same<T,float>::value, float  >::type s_b
                ,const RangeVector & y, float s_a,
                const CrsMatrix & A , const DomainVector & x) {

  char matdescra[6] = "GLNC0";
  int stride_x = layout_x.stride[0];
  int stride_y = layout_y.stride[0];
  char transa = 'N';
  int m = A.numRows();
  int n = x.dimension_1();
  int k = A.numCols();
  float* x_ptr = (float*)x.ptr_on_device();
  float* y_ptr = (float*)y.ptr_on_device();
  if(x.dimension_1()>1) {

    Impl::PhysicalLayout layout_x(x);
    Impl::PhysicalLayout layout_y(y);
    if((layout_x.layout_type!=layout_x.Right) || layout_y.layout_type!=layout_y.Right) return false;

    mkl_scsrmm(&transa,
               &m, &n, &k,
               &s_a,
               matdescra,
               A.values.ptr_on_device(),
               A.graph.entries.ptr_on_device(),
               (int*) &A.graph.row_map(0),
               (int*) &A.graph.row_map(1),
               x_ptr,
               &stride_x,
               &s_b,
               y_ptr,
               &stride_y);
  } else
    mkl_scsrmv(&transa,
              &m, &k,
              &s_a,
              matdescra,
              A.values.ptr_on_device(),
              A.graph.entries.ptr_on_device(),
              (int*) &A.graph.row_map(0),
              (int*) &A.graph.row_map(1),
              x_ptr,
              &s_b,
              y_ptr);
  return true;
}

//ToDo: strip compatible type attributes (const, volatile); make type of s_b and s_a independent
template<class RangeVector,class CrsMatrix,class DomainVector>
bool MV_Multiply_Try_MKL( typename RangeVector::value_type s_b,const RangeVector & y, typename DomainVector::value_type s_a,
                const CrsMatrix & A , const DomainVector & x)
{
  if( ! Kokkos::Impl::is_same<typename RangeVector::device_type::memory_space,typename Kokkos::HostSpace>::value ) return false;
  if(Kokkos::Impl::is_same<typename RangeVector::non_const_value_type,float>::value&&
         Kokkos::Impl::is_same<typename DomainVector::non_const_value_type,float>::value&&
         Kokkos::Impl::is_same<typename CrsMatrix::values_type::non_const_value_type,float>::value) {
           return MV_Multiply_DoMKL<typename RangeVector::value_type,RangeVector,CrsMatrix,DomainVector>(s_b,y,s_a,A,x);
  } else
  if(Kokkos::Impl::is_same<typename RangeVector::non_const_value_type,double>::value&&
         Kokkos::Impl::is_same<typename DomainVector::non_const_value_type,double>::value&&
         Kokkos::Impl::is_same<typename CrsMatrix::values_type::non_const_value_type,double>::value) {
           return MV_Multiply_DoMKL<typename RangeVector::value_type,RangeVector,CrsMatrix,DomainVector>(s_b,y,s_a,A,x);
  } else
  return false;
}

}

#endif // KOKKOS_CRSMATRIX_MKL_HPP
