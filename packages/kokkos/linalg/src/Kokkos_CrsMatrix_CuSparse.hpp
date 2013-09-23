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
#ifdef KOKKOS_CRSMATRIX_CUSPARSE_HPP_
#define KOKKOS_CRSMATRIX_CUSPARSE_HPP_

/// \file Kokkos_CrsMatrix_CuSparse.hpp
/// \brief cuSPARSE implementation of Kokkos::CrsMatrix kernels.

#include <impl/Kokkos_PhysicalLayout.hpp>

namespace Kokkos {
/// \namespace CuSparse
/// \brief cuSPARSE implementation of Kokkos kernels.
///
/// \warning This namespace and everything in it is a implementation
///   detail of Kokkos.  Do not rely on this namespace or anything
///   in it.  It may change or disappear at any time.
namespace CuSparse {

/// \fn MV_Multiply_DoCuSparse
/// \brief Attempt to use cuSPARSE to implement sparse matrix-(multi)vector multiply.
/// \tparam T The type of entries in the output (range) vector(s).
/// \tparam RangeVectorType The type of the output (range) vector(s).
/// \tparam CrsMatrixType The type of the sparse matrix; must be a
///   Kokkos::CrsMatrix specialization.
/// \tparam DomainVectorType The type of the input (domain) vector(s).
///
/// \return \c true if the attempt succeeded, else \c false.
template<typename T, class RangeVectorType, class CrsMatrixType, class DomainVectorType>
bool
MV_Multiply_DoCuSparse (typename Kokkos::Impl::enable_if<! Kokkos::Impl::is_same<T, double>::value && ! Kokkos::Impl::is_same<T, float>::value, typename RangeVectorType::scalar_type>::type s_b, const RangeVectorType& y, typename DomainVectorType::scalar_type s_a, const CrsMatrixType& A , const DomainVectorType& x) {
  return false;
}

template<typename T, class RangeVector,class CrsMatrix,class DomainVector>
bool
MV_Multiply_DoCuSparse (typename Kokkos::Impl::enable_if<Kokkos::Impl::is_same<T, double>::value, double>::type s_b,
                        const RangeVector& y,
                        double s_a,
                        const CrsMatrix& A,
                        const DomainVector& x)
{

  if (x.dimension_1 () == 1) {
    cusparseDcsrmv (A.cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                    A.numRows(), A.numCols(),  A.nnz(),
                    &s_a,
                    A.cusparse_descr,
                    A.values.ptr_on_device(),
                    (const int*) A.graph.row_map.ptr_on_device(),
                    A.graph.entries.ptr_on_device(),
                    x.ptr_on_device(),
                    &s_b,
                    y.ptr_on_device());
  } else {
    Impl::PhysicalLayout layout_x (x);
    Impl::PhysicalLayout layout_y (y);
    if ((layout_x.layout_type != layout_x.Left) || layout_y.layout_type != layout_y.Left) {
      return false;
    }
    cusparseDcsrmm (A.cusparse_handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
                    A.numRows(), x.dimension_1(), A.numCols(),  A.nnz(),
                    &s_a,
                    A.cusparse_descr,
                    A.values.ptr_on_device(),
                    (const int*) A.graph.row_map.ptr_on_device(),
                    A.graph.entries.ptr_on_device(),
                    x.ptr_on_device(),
                    layout_x.stride[1],
                    &s_b,
                    y.ptr_on_device(),
                    layout_y.stride[1]);
  }
  return true;
}

template<typename T, class RangeVector,class CrsMatrix,class DomainVector>
bool MV_Multiply_DoCuSparse(typename Kokkos::Impl::enable_if<Kokkos::Impl::is_same<T,float>::value, float  >::type s_b
                ,const RangeVector & y, float s_a,
                const CrsMatrix & A , const DomainVector & x) {
        if(x.dimension_1()==1) {
        cusparseScsrmv(A.cusparse_handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
                       A.numRows(), A.numCols(),  A.nnz(),
                       &s_a,
                       A.cusparse_descr,
                       A.values.ptr_on_device(),
                       (const int*) A.graph.row_map.ptr_on_device(),
                       A.graph.entries.ptr_on_device(),
                       x.ptr_on_device(),
                       &s_b,
                       y.ptr_on_device());
        } else {
          Impl::PhysicalLayout layout_x(x);
          Impl::PhysicalLayout layout_y(y);
          if((layout_x.layout_type!=layout_x.Left) || layout_y.layout_type!=layout_y.Left) return false;
          cusparseScsrmm(A.cusparse_handle,CUSPARSE_OPERATION_NON_TRANSPOSE,
                               A.numRows(), x.dimension_1(), A.numCols(),  A.nnz(),
                               &s_a,
                               A.cusparse_descr,
                               A.values.ptr_on_device(),
                               (const int*) A.graph.row_map.ptr_on_device(),
                               A.graph.entries.ptr_on_device(),
                               x.ptr_on_device(),
                               layout_x.stride[1],
                               &s_b,
                               y.ptr_on_device(),
                               layout_y.stride[1]);     
  }
  return true;
}

//ToDo: strip compatible type attributes (const, volatile); make type of s_b and s_a independent
template<class RangeVector,class CrsMatrix,class DomainVector>
bool
MV_Multiply_Try_CuSparse (typename RangeVector::scalar_type s_b,
                          const RangeVector& y,
                          typename DomainVector::scalar_type s_a,
                          const CrsMatrix& A,
                          const DomainVector& x)
{
  if(!Kokkos::Impl::is_same<typename RangeVector::device_type,typename Kokkos::Cuda>::value) return false;
  if(Kokkos::Impl::is_same<typename RangeVector::non_const_scalar_type,float>::value&&
         Kokkos::Impl::is_same<typename DomainVector::non_const_scalar_type,float>::value&&
         Kokkos::Impl::is_same<typename CrsMatrix::values_type::non_const_scalar_type,float>::value) {
           return MV_Multiply_DoCuSparse<typename RangeVector::scalar_type,RangeVector,CrsMatrix,DomainVector>(s_b,y,s_a,A,x);
  } else
  if(Kokkos::Impl::is_same<typename RangeVector::non_const_scalar_type,double>::value&&
         Kokkos::Impl::is_same<typename DomainVector::non_const_scalar_type,double>::value&&
         Kokkos::Impl::is_same<typename CrsMatrix::values_type::non_const_scalar_type,double>::value) {
           return MV_Multiply_DoCuSparse<typename RangeVector::scalar_type,RangeVector,CrsMatrix,DomainVector>(s_b,y,s_a,A,x);
  } else
  return false;
}

} // namespace CuSparse
} // namespace Kokkos

#endif // KOKKOS_CRSMATRIX_CUSPARSE_HPP_
