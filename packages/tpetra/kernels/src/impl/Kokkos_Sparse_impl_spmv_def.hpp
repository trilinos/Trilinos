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

#ifndef KOKKOS_SPARSE_IMPL_SPMV_DEF_HPP_
#define KOKKOS_SPARSE_IMPL_SPMV_DEF_HPP_

#include "Kokkos_InnerProductSpaceTraits.hpp"
#include "Kokkos_CrsMatrix.hpp" // only included for backwards compatibility
#include "Kokkos_Blas1_MV.hpp"

#ifdef HAVE_TPETRAKERNELS_ETI_ONLY
#define KOKKOSSPARSE_ETI_ONLY
#endif

namespace KokkosSparse {
namespace Impl {

// This TransposeFunctor is functional, but not necessarily performant.
template<class AMatrix,
         class XVector,
         class YVector,
         int doalpha,
         int dobeta,
         bool conjugate>
struct SPMV_Transpose_Functor {
  typedef typename AMatrix::execution_space            execution_space;
  typedef typename AMatrix::non_const_ordinal_type     ordinal_type;
  typedef typename AMatrix::non_const_value_type       value_type;
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type            team_member;
  typedef Kokkos::Details::ArithTraits<value_type>     ATV;

  const value_type alpha;
  AMatrix  m_A;
  XVector m_x;
  const value_type beta;
  YVector m_y;

  const ordinal_type rows_per_thread;

  SPMV_Transpose_Functor (const value_type alpha_,
                          const AMatrix m_A_,
                          const XVector m_x_,
                          const value_type beta_,
                          const YVector m_y_,
                          const ordinal_type rows_per_thread_) :
    alpha (alpha_), m_A (m_A_), m_x (m_x_),
    beta (beta_), m_y (m_y_),
    rows_per_thread (rows_per_thread_)
  {}

  KOKKOS_INLINE_FUNCTION void
  operator() (const team_member& dev) const
  {
    // This should be a thread loop as soon as we can use C++11
    for (ordinal_type loop = 0; loop < rows_per_thread; ++loop) {
      // iRow represents a row of the matrix, so its correct type is
      // ordinal_type.
      const ordinal_type iRow = (static_cast<ordinal_type> (dev.league_rank() * dev.team_size() + dev.team_rank()))
                                * rows_per_thread + loop;
      if (iRow >= m_A.numRows ()) {
        return;
      }

      const auto row = m_A.rowConst (iRow);
      const ordinal_type row_length = row.length;

#ifdef __CUDA_ARCH__
      for (ordinal_type iEntry = static_cast<ordinal_type> (threadIdx.x);
           iEntry < row_length;
           iEntry += static_cast<ordinal_type> (blockDim.x)) {
#else
      for (ordinal_type iEntry = 0;
           iEntry < row_length;
           iEntry ++) {
#endif
        const value_type val = conjugate ?
          ATV::conj (row.value(iEntry)) :
          row.value(iEntry);
        const ordinal_type ind = row.colidx(iEntry);

        if (doalpha != 1) {
          Kokkos::atomic_add (&m_y(ind), value_type(alpha * val * m_x(iRow)));
        } else {
          Kokkos::atomic_add (&m_y(ind), value_type(val * m_x(iRow)));
        }
      }
    }
  }
};

template<class AMatrix,
         class XVector,
         class YVector,
         int doalpha,
         int dobeta,
         bool conjugate>
struct SPMV_Functor {
  typedef typename AMatrix::execution_space            execution_space;
  typedef typename AMatrix::non_const_ordinal_type     ordinal_type;
  typedef typename AMatrix::non_const_value_type       value_type;
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type            team_member;
  typedef Kokkos::Details::ArithTraits<value_type>     ATV;

  const value_type alpha;
  AMatrix  m_A;
  XVector m_x;
  const value_type beta;
  YVector m_y;

  const ordinal_type rows_per_thread;

  SPMV_Functor (const value_type alpha_,
               const AMatrix m_A_,
               const XVector m_x_,
               const value_type beta_,
               const YVector m_y_,
               const ordinal_type rows_per_thread_) :
    alpha (alpha_), m_A (m_A_), m_x (m_x_),
    beta (beta_), m_y (m_y_),
    rows_per_thread (rows_per_thread_)
  {}

  KOKKOS_INLINE_FUNCTION void
  operator() (const team_member& dev) const
  {
    // This should be a thread loop as soon as we can use C++11
    for (ordinal_type loop = 0; loop < rows_per_thread; ++loop) {
      // iRow represents a row of the matrix, so its correct type is
      // ordinal_type.
      const ordinal_type iRow = (static_cast<ordinal_type> (dev.league_rank() * dev.team_size() + dev.team_rank()))
                                * rows_per_thread + loop;
      if (iRow >= m_A.numRows ()) {
        return;
      }
      const auto row = m_A.rowConst (iRow);
      const ordinal_type row_length = row.length;
      value_type sum = 0;

      // Use explicit Cuda below to avoid C++11 for now. This should be a vector reduce loop !
      #ifdef KOKKOS_HAVE_PRAGMA_IVDEP
      #pragma ivdep
      #endif
      #ifdef KOKKOS_HAVE_PRAGMA_UNROLL
      #pragma unroll
      #endif
      #ifdef KOKKOS_HAVE_PRAGMA_LOOPCOUNT
      #pragma loop count (15)
      #endif
#ifdef __CUDA_ARCH__
      for (ordinal_type iEntry = static_cast<ordinal_type> (threadIdx.x);
           iEntry < row_length;
           iEntry += static_cast<ordinal_type> (blockDim.x)) {
#else
      for (ordinal_type iEntry = 0;
           iEntry < row_length;
           iEntry ++) {
#endif
        const value_type val = conjugate ?
                ATV::conj (row.value(iEntry)) :
                row.value(iEntry);
        sum += val * m_x(row.colidx(iEntry));
      }

#if defined(__CUDA_ARCH__) && defined(KOKKOS_HAVE_CUDA)
      if (blockDim.x > 1)
        sum += Kokkos::shfl_down(sum, 1,blockDim.x);
      if (blockDim.x > 2)
        sum += Kokkos::shfl_down(sum, 2,blockDim.x);
      if (blockDim.x > 4)
        sum += Kokkos::shfl_down(sum, 4,blockDim.x);
      if (blockDim.x > 8)
        sum += Kokkos::shfl_down(sum, 8,blockDim.x);
      if (blockDim.x > 16)
        sum += Kokkos::shfl_down(sum, 16,blockDim.x);

      if (threadIdx.x==0) {
#else
      if (true) {
#endif // defined(__CUDA_ARCH__) && defined(KOKKOS_HAVE_CUDA)
        if (doalpha == -1) {
          sum *= value_type(-1);
        } else if (doalpha * doalpha != 1) {
          sum *= alpha;
        }

        if (dobeta == 0) {
          m_y(iRow) = sum ;
        } else if (dobeta == 1) {
          m_y(iRow) += sum ;
        } else if (dobeta == -1) {
          m_y(iRow) = -m_y(iRow) +  sum;
        } else {
          m_y(iRow) = beta * m_y(iRow) + sum;
        }
      }
    }
  }
};

template<class AMatrix,
         class XVector,
         class YVector,
         int doalpha,
         int dobeta,
         bool conjugate>
static void spmv_alpha_beta_no_transpose(typename AMatrix::const_value_type& alpha, const AMatrix& A, const XVector& x, typename AMatrix::const_value_type& beta, const YVector& y) {
  typedef typename AMatrix::ordinal_type ordinal_type;

  if (A.numRows () <= static_cast<ordinal_type> (0)) {
    return;
  }
  if (doalpha == 0) {
    if (dobeta==2) {
      KokkosBlas::scal(y, beta, y);
    }
    else {
      KokkosBlas::scal(y, typename YVector::const_value_type(beta), y);
    }
    return;
  } else {
    typedef typename AMatrix::size_type size_type;

    //Impl::MV_Multiply_Check_Compatibility(betav,y,alphav,A,x,doalpha,dobeta);

    // Assuming that no row contains duplicate entries, NNZPerRow
    // cannot be more than the number of columns of the matrix.  Thus,
    // the appropriate type is ordinal_type.
    const ordinal_type NNZPerRow = static_cast<ordinal_type> (A.nnz () / A.numRows ());

    int vector_length = 1;
    while( (static_cast<ordinal_type> (vector_length*2*3) <= NNZPerRow) && (vector_length<32) ) vector_length*=2;

#ifndef KOKKOS_FAST_COMPILE // This uses templated functions on doalpha and dobeta and will produce 16

    typedef SPMV_Functor<AMatrix, XVector, YVector, doalpha, dobeta , conjugate> OpType ;

    typename AMatrix::const_ordinal_type nrow = A.numRows();

    OpType op(alpha,A,x,beta,y,RowsPerThread<typename AMatrix::execution_space >(NNZPerRow)) ;

    const int rows_per_thread = RowsPerThread<typename AMatrix::execution_space >(NNZPerRow);
    const int team_size = Kokkos::TeamPolicy< typename AMatrix::execution_space >::team_size_recommended(op,vector_length);
    const int rows_per_team = rows_per_thread * team_size;
    const size_type nteams = (nrow+rows_per_team-1)/rows_per_team;
    Kokkos::parallel_for( Kokkos::TeamPolicy< typename AMatrix::execution_space , Kokkos::Schedule<Kokkos::Dynamic> >
       ( nteams , team_size , vector_length ) , op );

#else // KOKKOS_FAST_COMPILE this will only instantiate one Kernel for alpha/beta

    typedef SPMV_Functor<AMatrix, XVector, YVector, 2, 2 , conjugate> OpType ;

    typename AMatrix::const_ordinal_type nrow = A.numRows();

    OpType op(alpha,A,x,beta,y,RowsPerThread<typename AMatrix::execution_space >(NNZPerRow)) ;

    const int rows_per_thread = RowsPerThread<typename AMatrix::execution_space >(NNZPerRow);
    const int team_size = Kokkos::TeamPolicy< typename AMatrix::execution_space >::team_size_recommended(op,vector_length);
    const int rows_per_team = rows_per_thread * team_size;
    const size_type nteams = (nrow+rows_per_team-1)/rows_per_team;
    Kokkos::parallel_for( Kokkos::TeamPolicy< typename AMatrix::execution_space >
       ( nteams , team_size , vector_length ) , op );

#endif // KOKKOS_FAST_COMPILE
  }
}

template<class AMatrix,
         class XVector,
         class YVector,
         int doalpha,
         int dobeta,
         bool conjugate>
static void spmv_alpha_beta_transpose(typename AMatrix::const_value_type& alpha, const AMatrix& A, const XVector& x, typename AMatrix::const_value_type& beta, const YVector& y) {
  typedef typename AMatrix::ordinal_type ordinal_type;

  if (A.numRows () <= static_cast<ordinal_type> (0)) {
    return;
  }

  // Scale y since we atomic add to it.
  if (dobeta==2) {
    KokkosBlas::scal(y, beta, y);
  }
  else {
    KokkosBlas::scal(y, typename YVector::const_value_type(beta), y);
  }

  if (doalpha != 0) {
    typedef typename AMatrix::size_type size_type;

    //Impl::MV_Multiply_Check_Compatibility(betav,y,alphav,A,x,doalpha,dobeta);

    // Assuming that no row contains duplicate entries, NNZPerRow
    // cannot be more than the number of columns of the matrix.  Thus,
    // the appropriate type is ordinal_type.
    const ordinal_type NNZPerRow = static_cast<ordinal_type> (A.nnz () / A.numRows ());

    int vector_length = 1;
    while( (static_cast<ordinal_type> (vector_length*2*3) <= NNZPerRow) && (vector_length<32) ) vector_length*=2;

#ifndef KOKKOS_FAST_COMPILE // This uses templated functions on doalpha and dobeta and will produce 16

    typedef SPMV_Transpose_Functor<AMatrix, XVector, YVector, doalpha, dobeta , conjugate> OpType ;

    typename AMatrix::const_ordinal_type nrow = A.numRows();

    OpType op(alpha,A,x,beta,y,RowsPerThread<typename AMatrix::execution_space >(NNZPerRow)) ;

    const int rows_per_thread = RowsPerThread<typename AMatrix::execution_space >(NNZPerRow);
    const int team_size = Kokkos::TeamPolicy< typename AMatrix::execution_space >::team_size_recommended(op,vector_length);
    const int rows_per_team = rows_per_thread * team_size;
    const size_type nteams = (nrow+rows_per_team-1)/rows_per_team;
    Kokkos::parallel_for( Kokkos::TeamPolicy< typename AMatrix::execution_space >
       ( nteams , team_size , vector_length ) , op );

#else // KOKKOS_FAST_COMPILE this will only instantiate one Kernel for alpha/beta

    typedef SPMV_Transpose_Functor<AMatrix, XVector, YVector, 2, 2 , conjugate> OpType ;

    typename AMatrix::const_ordinal_type nrow = A.numRows();

    OpType op(alpha,A,x,beta,y,RowsPerThread<typename AMatrix::execution_space >(NNZPerRow)) ;

    const int rows_per_thread = RowsPerThread<typename AMatrix::execution_space >(NNZPerRow);
    const int team_size = Kokkos::TeamPolicy< typename AMatrix::execution_space >::team_size_recommended(op,vector_length);
    const int rows_per_team = rows_per_thread * team_size;
    const size_type nteams = (nrow+rows_per_team-1)/rows_per_team;
    Kokkos::parallel_for( Kokkos::TeamPolicy< typename AMatrix::execution_space >
       ( nteams , team_size , vector_length ) , op );

#endif // KOKKOS_FAST_COMPILE
  }
}

template<class AMatrix,
         class XVector,
         class YVector,
         int doalpha,
         int dobeta>
static void spmv_alpha_beta(const char mode[], typename AMatrix::const_value_type& alpha, const AMatrix& A, const XVector& x, typename AMatrix::const_value_type& beta, const YVector& y) {
  if(mode[0]==NoTranspose[0]) {
    spmv_alpha_beta_no_transpose<AMatrix,XVector,YVector,doalpha,dobeta,false>
      (alpha,A,x,beta,y);
    return;
  }
  if(mode[0]==Conjugate[0]) {
    spmv_alpha_beta_no_transpose<AMatrix,XVector,YVector,doalpha,dobeta,true>
      (alpha,A,x,beta,y);
    return;
  }
  if(mode[0]==Transpose[0]) {
    spmv_alpha_beta_transpose<AMatrix,XVector,YVector,doalpha,dobeta,false>
      (alpha,A,x,beta,y);
    return;
  }
  if(mode[0]==ConjugateTranspose[0]) {
    spmv_alpha_beta_transpose<AMatrix,XVector,YVector,doalpha,dobeta,true>
      (alpha,A,x,beta,y);
    return;
  }
  Kokkos::Impl::throw_runtime_exception("Invalid Transpose Mode for KokkosSparse::spmv()");
}

template<class AMatrix,
         class XVector,
         class YVector,
         int doalpha>
void spmv_alpha(const char mode[], typename AMatrix::const_value_type& alpha, const AMatrix& A, const XVector& x, typename AMatrix::const_value_type& beta, const YVector& y) {
  typedef typename AMatrix::non_const_value_type Scalar;

  if( beta == Kokkos::Details::ArithTraits<Scalar>::zero () ) {
    spmv_alpha_beta<AMatrix,XVector,YVector,doalpha,0>(mode,alpha,A,x,beta,y);
    return;
  }
  if( beta == Kokkos::Details::ArithTraits<Scalar>::one () ) {
    spmv_alpha_beta<AMatrix,XVector,YVector,doalpha,1>(mode,alpha,A,x,beta,y);
    return;
  }
  if( beta == -Kokkos::Details::ArithTraits<Scalar>::one () ) {
    spmv_alpha_beta<AMatrix,XVector,YVector,doalpha,-1>(mode,alpha,A,x,beta,y);
    return;
  }
  spmv_alpha_beta<AMatrix,XVector,YVector,doalpha,2>(mode,alpha,A,x,beta,y);
}

#ifndef KOKKOSSPARSE_ETI_ONLY
template<class AT, class AO, class AD, class AM, class AS,
         class XT, class XL, class XD, class XM,
         class YT, class YL, class YD, class YM>
void
SPMV<AT, AO, AD, AM, AS,
     XT, XL, XD, XM,
     YT, YL, YD, YM>::
spmv (const char mode[], const Scalar& alpha, const AMatrix& A,
      const XVector& x, const Scalar& beta, const YVector& y)
{
  if (alpha == Kokkos::Details::ArithTraits<Scalar>::zero ()) {
    spmv_alpha<AMatrix,XVector,YVector,0> (mode, alpha, A, x, beta, y);
    return;
  }
  if (alpha == Kokkos::Details::ArithTraits<Scalar>::one ()) {
    spmv_alpha<AMatrix,XVector,YVector,1> (mode, alpha, A, x, beta, y);
    return;
  }
  if (alpha == -Kokkos::Details::ArithTraits<Scalar>::one ()) {
    spmv_alpha<AMatrix,XVector,YVector,-1> (mode, alpha, A, x, beta, y);
    return;
  }
  spmv_alpha<AMatrix,XVector,YVector,2> (mode, alpha, A, x, beta, y);
}
#endif // KOKKOSSPARSE_ETI_ONLY

//
// Macro for defining (not declaring; for the declaration macro, see
// above) a full specialization of the SPMV struct.  We use this macro
// in various .cpp file(s) in this directory.
//
// SCALAR_TYPE: The type of each entry in the sparse matrix
// ORDINAL_TYPE: The type of each column index in the sparse matrix
// OFFSET_TYPE: The type of each row offset in the sparse matrix
// LAYOUT_TYPE: The layout of the Kokkos::View vector arguments
//   of the sparse matrix-vector multiply
// EXEC_SPACE_TYPE: The execution space type
// MEM_SPACE_TYPE: The memory space type
//
#define KOKKOSSPARSE_IMPL_SPMV_DEF( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE ) \
void \
SPMV<const SCALAR_TYPE, \
     ORDINAL_TYPE, \
     Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
     Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
     OFFSET_TYPE, \
     const SCALAR_TYPE*, \
     Kokkos::LayoutLeft, \
     Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
     Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess>, \
     SCALAR_TYPE*, \
     Kokkos::LayoutLeft, \
     Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
     Kokkos::MemoryTraits<Kokkos::Unmanaged> >:: \
spmv (const char mode[], const Scalar& alpha, const AMatrix& A, \
      const XVector& x, const Scalar& beta, const YVector& y) \
{ \
  if (alpha == Kokkos::Details::ArithTraits<Scalar>::zero ()) { \
    spmv_alpha<AMatrix,XVector,YVector,0> (mode, alpha, A, x, beta, y); \
    return; \
  } \
  if (alpha == Kokkos::Details::ArithTraits<Scalar>::one ()) { \
    spmv_alpha<AMatrix,XVector,YVector,1> (mode, alpha, A, x, beta, y); \
    return; \
  } \
  if (alpha == -Kokkos::Details::ArithTraits<Scalar>::one ()) { \
    spmv_alpha<AMatrix,XVector,YVector,-1> (mode, alpha, A, x, beta, y); \
    return; \
  } \
  spmv_alpha<AMatrix,XVector,YVector,2> (mode, alpha, A, x, beta, y); \
}

//
// Macro for defining (not declaring; for the declaration macro, see
// above) a full specialization of the SPMV struct, using the default
// OFFSET_TYPE and LAYOUT_TYPE = Kokkos::LayoutLeft (Tpetra's layout).
// We use this macro in various .cpp file(s) in this directory.
//
// SCALAR_TYPE: The type of each entry in the sparse matrix
// ORDINAL_TYPE: The type of each column index in the sparse matrix
// EXEC_SPACE_TYPE: The execution space type
// MEM_SPACE_TYPE: The memory space type
//
#define KOKKOSSPARSE_IMPL_SPMV_DEFAULTS_DEF( SCALAR_TYPE, ORDINAL_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE ) \
void \
SPMV<const SCALAR_TYPE, \
     ORDINAL_TYPE, \
     Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
     Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
     Kokkos::StaticCrsGraph<ORDINAL_TYPE, Kokkos::LayoutLeft, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE> >::size_type, \
     const SCALAR_TYPE*, \
     Kokkos::LayoutLeft, \
     Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
     Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess>, \
     SCALAR_TYPE*, \
     Kokkos::LayoutLeft, \
     Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
     Kokkos::MemoryTraits<Kokkos::Unmanaged> >:: \
spmv (const char mode[], const Scalar& alpha, const AMatrix& A, \
      const XVector& x, const Scalar& beta, const YVector& y) \
{ \
  if (alpha == Kokkos::Details::ArithTraits<Scalar>::zero ()) { \
    spmv_alpha<AMatrix,XVector,YVector,0> (mode, alpha, A, x, beta, y); \
    return; \
  } \
  if (alpha == Kokkos::Details::ArithTraits<Scalar>::one ()) { \
    spmv_alpha<AMatrix,XVector,YVector,1> (mode, alpha, A, x, beta, y); \
    return; \
  } \
  if (alpha == -Kokkos::Details::ArithTraits<Scalar>::one ()) { \
    spmv_alpha<AMatrix,XVector,YVector,-1> (mode, alpha, A, x, beta, y); \
    return; \
  } \
  spmv_alpha<AMatrix,XVector,YVector,2> (mode, alpha, A, x, beta, y); \
}


// Functor for implementing transpose and conjugate transpose sparse
// matrix-vector multiply with multivector (2-D View) input and
// output.  This functor works, but is not necessarily performant.
template<class aCoeffs,
         class AMatrix,
         class XVector,
         class bCoeffs,
         class YVector,
         int doalpha,
         int dobeta,
         bool conjugate>
struct SPMV_MV_Transpose_Functor {
  typedef typename AMatrix::execution_space            execution_space;
  typedef typename AMatrix::non_const_ordinal_type     ordinal_type;
  typedef typename AMatrix::non_const_value_type       value_type;
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type            team_member;
  typedef Kokkos::Details::ArithTraits<value_type>     ATV;

  aCoeffs alpha;
  AMatrix  m_A;
  XVector  m_x;
  bCoeffs beta;
  YVector  m_y;

  const ordinal_type n;
  const ordinal_type rows_per_thread;

  SPMV_MV_Transpose_Functor (const aCoeffs alpha_,
                             const AMatrix m_A_,
                             const XVector m_x_,
                             const bCoeffs beta_,
                             const YVector m_y_,
                             const ordinal_type rows_per_thread_) :
    alpha (alpha_),
    m_A (m_A_), m_x (m_x_), beta (beta_), m_y (m_y_), n (m_x_.dimension_1()),
    rows_per_thread (rows_per_thread_)
  {}

  KOKKOS_INLINE_FUNCTION void
  operator() (const team_member& dev) const
  {
    // This should be a thread loop as soon as we can use C++11
    for (ordinal_type loop = 0; loop < rows_per_thread; ++loop) {
      // iRow represents a row of the matrix, so its correct type is
      // ordinal_type.
      const ordinal_type iRow = (static_cast<ordinal_type> (dev.league_rank() * dev.team_size() + dev.team_rank()))
                                * rows_per_thread + loop;
      if (iRow >= m_A.numRows ()) {
        return;
      }

      const auto row = m_A.rowConst (iRow);
      const ordinal_type row_length = row.length;

#ifdef __CUDA_ARCH__
      for (ordinal_type iEntry = static_cast<ordinal_type> (threadIdx.x);
           iEntry < static_cast<ordinal_type> (row_length);
           iEntry += static_cast<ordinal_type> (blockDim.x)) {
#else
      for (ordinal_type iEntry = 0;
           iEntry < row_length;
           iEntry ++) {
#endif
        const value_type val = conjugate ?
          ATV::conj (row.value(iEntry)) :
          row.value(iEntry);
        const ordinal_type ind = row.colidx(iEntry);

        if (doalpha != 1) {
          #ifdef KOKKOS_HAVE_PRAGMA_UNROLL
          #pragma unroll
          #endif
          for (ordinal_type k = 0; k < n; ++k) {
            Kokkos::atomic_add (&m_y(ind,k), value_type(alpha(k) * val * m_x(iRow, k)));
          }
        } else {
          #ifdef KOKKOS_HAVE_PRAGMA_UNROLL
          #pragma unroll
          #endif
          for (ordinal_type k = 0; k < n; ++k) {
            Kokkos::atomic_add (&m_y(ind,k), value_type(val * m_x(iRow, k)));
          }
        }
      }
    }
  }
};

template<class aCoeffs,
         class AMatrix,
         class XVector,
         class bCoeffs,
         class YVector,
         int doalpha,
         int dobeta,
         bool conjugate>
struct SPMV_MV_LayoutLeft_Functor {
  typedef typename AMatrix::execution_space      execution_space;
  typedef typename AMatrix::ordinal_type         ordinal_type;
  typedef typename AMatrix::non_const_value_type value_type;
  typedef typename Kokkos::TeamPolicy<execution_space>        team_policy;
  typedef typename team_policy::member_type                   team_member;
  typedef Kokkos::Details::ArithTraits<value_type>     ATV;

  aCoeffs alpha;
  AMatrix  m_A;
  XVector  m_x;
  bCoeffs beta;
  YVector  m_y;
  //! The number of columns in the input and output MultiVectors.
  ordinal_type n;
  ordinal_type rows_per_thread;

  SPMV_MV_LayoutLeft_Functor (const aCoeffs alpha_,
                              const AMatrix m_A_,
                              const XVector m_x_,
                              const bCoeffs beta_,
                              const YVector m_y_,
                              const ordinal_type rows_per_thread_) :
    alpha (alpha_),
    m_A (m_A_), m_x (m_x_), beta (beta_), m_y (m_y_), n (m_x_.dimension_1()),
    rows_per_thread (rows_per_thread_)
  {}

  template<int UNROLL>
  KOKKOS_INLINE_FUNCTION void
  strip_mine (const team_member& dev, const ordinal_type& iRow, const ordinal_type& kk) const
  {
    value_type sum[UNROLL];

#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
    for (int k = 0; k < UNROLL; ++k) {
      // NOTE (mfh 09 Aug 2013) This requires that assignment from int
      // (in this case, 0) to value_type be defined.  It's not for
      // types like arprec and dd_real.
      //
      // mfh 29 Sep 2013: On the other hand, arprec and dd_real won't
      // work on CUDA devices anyway, since their methods aren't
      // device functions.  arprec has other issues (e.g., dynamic
      // memory allocation, and the array-of-structs memory layout
      // which is unfavorable to GPUs), but could be dealt with in the
      // same way as Sacado's AD types.
      sum[k] = 0;
    }

    const auto row = m_A.rowConst (iRow);

    // NOTE (mfh 20 Mar 2015) Unfortunately, Kokkos::Vectorization
    // lacks a typedef for determining the type of the return value of
    // begin().  I know that it returns int now, but this may change
    // at some point.
    //
    // The correct type of iEntry is ordinal_type.  This is because we
    // assume that rows have no duplicate entries.  As a result, a row
    // cannot have more entries than the number of columns in the
    // matrix.

#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
#ifdef KOKKOS_HAVE_PRAGMA_LOOPCOUNT
#pragma loop count (15)
#endif
#ifdef __CUDA_ARCH__
        for (ordinal_type iEntry = static_cast<ordinal_type> (threadIdx.x);
             iEntry < row.length;
             iEntry += static_cast<ordinal_type> (blockDim.x)) {
#else
        for (ordinal_type iEntry = 0;
             iEntry < row.length;
             iEntry ++) {
#endif
      const value_type val = conjugate ?
                  ATV::conj (row.value(iEntry)) :
                  row.value(iEntry);
      const ordinal_type ind = row.colidx(iEntry);

#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        sum[k] +=  val * m_x(ind, kk + k);
      }
    }

    if (doalpha == -1) {
      for (int ii=0; ii < UNROLL; ++ii) {
        value_type sumt=sum[ii];
#if defined(__CUDA_ARCH__) && defined(KOKKOS_HAVE_CUDA)
        if (blockDim.x > 1)
          sumt += Kokkos::shfl_down(sumt, 1,blockDim.x);
        if (blockDim.x > 2)
          sumt += Kokkos::shfl_down(sumt, 2,blockDim.x);
        if (blockDim.x > 4)
          sumt += Kokkos::shfl_down(sumt, 4,blockDim.x);
        if (blockDim.x > 8)
          sumt += Kokkos::shfl_down(sumt, 8,blockDim.x);
        if (blockDim.x > 16)
          sumt += Kokkos::shfl_down(sumt, 16,blockDim.x);
#endif // defined(__CUDA_ARCH__) && defined(KOKKOS_HAVE_CUDA)
        sum[ii] = - sumt;
      }
    }
    else {
      for (int ii=0; ii < UNROLL; ++ii) {
        value_type sumt = sum[ii];
#if defined(__CUDA_ARCH__) && defined(KOKKOS_HAVE_CUDA)
        if (blockDim.x > 1)
          sumt += Kokkos::shfl_down(sumt, 1,blockDim.x);
        if (blockDim.x > 2)
          sumt += Kokkos::shfl_down(sumt, 2,blockDim.x);
        if (blockDim.x > 4)
          sumt += Kokkos::shfl_down(sumt, 4,blockDim.x);
        if (blockDim.x > 8)
          sumt += Kokkos::shfl_down(sumt, 8,blockDim.x);
        if (blockDim.x > 16)
          sumt += Kokkos::shfl_down(sumt, 16,blockDim.x);
#endif // defined(__CUDA_ARCH__) && defined(KOKKOS_HAVE_CUDA)
        sum[ii] = sumt;
      }
    }

#if defined(__CUDA_ARCH__) && defined(KOKKOS_HAVE_CUDA)
    if (threadIdx.x==0) {
#else
    if (true) {
#endif // defined(__CUDA_ARCH__) && defined(KOKKOS_HAVE_CUDA)
      if (doalpha * doalpha != 1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          sum[k] *= alpha(kk + k);
        }
      }

      if (dobeta == 0) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(iRow, kk + k) = sum[k];
        }
      } else if (dobeta == 1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(iRow, kk + k) += sum[k];
        }
      } else if (dobeta == -1) {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(iRow, kk + k) = -m_y(iRow, kk + k) +  sum[k];
        }
      } else {
#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(iRow, kk + k) = beta(kk + k) * m_y(iRow, kk + k) + sum[k] ;
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION void
  strip_mine_1 (const team_member& dev, const ordinal_type& iRow) const
  {
    value_type sum = 0;

    const auto row = m_A.rowConst (iRow);

    // NOTE (mfh 20 Mar 2015) Unfortunately, Kokkos::Vectorization
    // lacks a typedef for determining the type of the return value of
    // begin().  I know that it returns int now, but this may change
    // at some point.
    //
    // The correct type of iEntry is ordinal_type.  This is because we
    // assume that rows have no duplicate entries.  As a result, a row
    // cannot have more entries than the number of columns in the
    // matrix.

#ifdef KOKKOS_HAVE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_HAVE_PRAGMA_UNROLL
#pragma unroll
#endif
#ifdef KOKKOS_HAVE_PRAGMA_LOOPCOUNT
#pragma loop count (15)
#endif
#ifdef __CUDA_ARCH__
    for (ordinal_type iEntry = static_cast<ordinal_type> (threadIdx.x);
         iEntry < row.length;
         iEntry += static_cast<ordinal_type> (blockDim.x)) {
#else
    for (ordinal_type iEntry = 0;
         iEntry < row.length;
         iEntry ++) {
#endif
      const value_type val = conjugate ?
                  ATV::conj (row.value(iEntry)) :
                  row.value(iEntry);
      sum += val * m_x(row.colidx(iEntry),0);
    }
#if defined(__CUDA_ARCH__) && defined(KOKKOS_HAVE_CUDA)
    if (blockDim.x > 1)
      sum += Kokkos::shfl_down(sum, 1,blockDim.x);
    if (blockDim.x > 2)
      sum += Kokkos::shfl_down(sum, 2,blockDim.x);
    if (blockDim.x > 4)
      sum += Kokkos::shfl_down(sum, 4,blockDim.x);
    if (blockDim.x > 8)
      sum += Kokkos::shfl_down(sum, 8,blockDim.x);
    if (blockDim.x > 16)
      sum += Kokkos::shfl_down(sum, 16,blockDim.x);
#endif // defined(__CUDA_ARCH__) && defined(KOKKOS_HAVE_CUDA)

#if defined(__CUDA_ARCH__) && defined(KOKKOS_HAVE_CUDA)
    if (threadIdx.x==0) {
#else
    if (true) {
#endif // defined(__CUDA_ARCH__) && defined(KOKKOS_HAVE_CUDA)
      if (doalpha == -1) {
        sum *= value_type(-1);
      } else if (doalpha * doalpha != 1) {
        sum *= alpha(0);
      }

      if (dobeta == 0) {
        m_y(iRow, 0) = sum ;
      } else if (dobeta == 1) {
        m_y(iRow, 0) += sum ;
      } else if (dobeta == -1) {
        m_y(iRow, 0) = -m_y(iRow, 0) +  sum;
      } else {
        m_y(iRow, 0) = beta(0) * m_y(iRow, 0) + sum;
      }
    }
  }


  KOKKOS_INLINE_FUNCTION void
  operator() (const team_member& dev) const
  {
    for (ordinal_type loop = 0; loop < rows_per_thread; ++loop) {

      // NOTE (mfh 20 Mar 2015) Unfortunately, Kokkos::Vectorization
      // lacks a typedef for determining the type of the return value
      // of global_thread_rank().  I know that it returns int now, but
      // this may change at some point.
      //
      // iRow represents a row of the matrix, so its correct type is
      // ordinal_type.

      const ordinal_type iRow = (dev.league_rank() * dev.team_size() + dev.team_rank())
                                * rows_per_thread + loop;
      if (iRow >= m_A.numRows ()) {
        return;
      }

      // mfh 20 Mar 2015, 07 Jun 2016: This is ordinal_type because it
      // needs to have the same type as n.
      ordinal_type kk = 0;

#ifdef KOKKOS_FAST_COMPILE
      for (; kk + 4 <= n; kk += 4) {
        strip_mine<4>(dev, iRow, kk);
      }
      for( ; kk < n; ++kk) {
        strip_mine<1>(dev, iRow, kk);
      }
#else
#  ifdef __CUDA_ARCH__
      if ((n > 8) && (n % 8 == 1)) {
        strip_mine<9>(dev, iRow, kk);
        kk += 9;
      }
      for(; kk + 8 <= n; kk += 8)
        strip_mine<8>(dev, iRow, kk);
      if(kk < n)
        switch(n - kk) {
#  else // NOT a CUDA device
          if ((n > 16) && (n % 16 == 1)) {
            strip_mine<17>(dev, iRow, kk);
            kk += 17;
          }

          for (; kk + 16 <= n; kk += 16) {
            strip_mine<16>(dev, iRow, kk);
          }

          if(kk < n)
            switch(n - kk) {
            case 15:
              strip_mine<15>(dev, iRow, kk);
              break;

            case 14:
              strip_mine<14>(dev, iRow, kk);
              break;

            case 13:
              strip_mine<13>(dev, iRow, kk);
              break;

            case 12:
              strip_mine<12>(dev, iRow, kk);
              break;

            case 11:
              strip_mine<11>(dev, iRow, kk);
              break;

            case 10:
              strip_mine<10>(dev, iRow, kk);
              break;

            case 9:
              strip_mine<9>(dev, iRow, kk);
              break;

            case 8:
              strip_mine<8>(dev, iRow, kk);
              break;
#  endif // __CUDA_ARCH__
            case 7:
              strip_mine<7>(dev, iRow, kk);
              break;

            case 6:
              strip_mine<6>(dev, iRow, kk);
              break;

            case 5:
              strip_mine<5>(dev, iRow, kk);
              break;

            case 4:
              strip_mine<4>(dev, iRow, kk);
              break;

            case 3:
              strip_mine<3>(dev, iRow, kk);
              break;

            case 2:
              strip_mine<2>(dev, iRow, kk);
              break;

            case 1:
              strip_mine_1(dev, iRow);
              break;
            }
#endif // KOKKOS_FAST_COMPILE
        }
    }
  };


template<class aCoeffs,
         class AMatrix,
         class XVector,
         class bCoeffs,
         class YVector,
         int doalpha,
         int dobeta,
         bool conjugate>
static void spmv_alpha_beta_mv_no_transpose(aCoeffs alpha, const AMatrix& A, const XVector& x, bCoeffs beta, const YVector& y) {
  typedef typename AMatrix::ordinal_type ordinal_type;

  if (A.numRows () <= static_cast<ordinal_type> (0)) {
    return;
  }
  if (doalpha == 0) {
    if (dobeta==2) {
      KokkosBlas::scal(y, beta, y);
    }
    else {
      KokkosBlas::scal(y, typename YVector::value_type (dobeta), y);
    }
    return;
  } else {
    typedef typename AMatrix::size_type size_type;

    //Impl::MV_Multiply_Check_Compatibility(betav,y,alphav,A,x,doalpha,dobeta);

    // Assuming that no row contains duplicate entries, NNZPerRow
    // cannot be more than the number of columns of the matrix.  Thus,
    // the appropriate type is ordinal_type.
    const ordinal_type NNZPerRow = static_cast<ordinal_type> (A.nnz () / A.numRows ());

    int vector_length = 1;
    while( (static_cast<ordinal_type> (vector_length*2*3) <= NNZPerRow) && (vector_length<8) ) vector_length*=2;

#ifndef KOKKOS_FAST_COMPILE // This uses templated functions on doalpha and dobeta and will produce 16 kernels

    typedef SPMV_MV_LayoutLeft_Functor<aCoeffs,AMatrix, XVector,
                                       bCoeffs, YVector,
                                       doalpha, dobeta, conjugate > OpType ;
    OpType op(alpha,A,x,beta,y,RowsPerThread<typename AMatrix::execution_space >(NNZPerRow)) ;

    typename AMatrix::const_ordinal_type nrow = A.numRows();

    // FIXME (mfh 07 Jun 2016) Shouldn't we use ordinal_type here
    // instead of int?  For example, if the number of threads is 1,
    // then this is just the number of rows.  Ditto for rows_per_team.
    // team_size is a hardware resource thing so it might legitimately
    // be int.
    const int rows_per_thread = RowsPerThread<typename AMatrix::execution_space >(NNZPerRow);
    const int team_size = Kokkos::TeamPolicy< typename AMatrix::execution_space >::team_size_recommended(op,vector_length);
    const int rows_per_team = rows_per_thread * team_size;
    const size_type nteams = (nrow+rows_per_team-1)/rows_per_team;
    Kokkos::parallel_for( Kokkos::TeamPolicy< typename AMatrix::execution_space >
       ( nteams , team_size , vector_length ) , op );

#else // KOKKOS_FAST_COMPILE this will only instantiate one Kernel for alpha/beta

    {
      typedef typename aCoeffs::non_const_value_type Scalar;
      if(doalpha==0)
        alpha = GetCoeffView<Scalar,typename aCoeffs::device_type>::get_view(
            Kokkos::Details::ArithTraits<Scalar>::zero(),x.dimension_1());
      if(doalpha==1)
        alpha = GetCoeffView<Scalar,typename aCoeffs::device_type>::get_view(
            Kokkos::Details::ArithTraits<Scalar>::one(),x.dimension_1());
      if(doalpha==1)
        alpha = GetCoeffView<Scalar,typename aCoeffs::device_type>::get_view(
            -Kokkos::Details::ArithTraits<Scalar>::one(),x.dimension_1());
    }

    {
      typedef typename bCoeffs::non_const_value_type Scalar;
      if(dobeta==0)
        beta = GetCoeffView<Scalar,typename aCoeffs::device_type>::get_view(
            Kokkos::Details::ArithTraits<Scalar>::zero(),x.dimension_1());
      if(dobeta==1)
        beta = GetCoeffView<Scalar,typename aCoeffs::device_type>::get_view(
            Kokkos::Details::ArithTraits<Scalar>::one(),x.dimension_1());
      if(dobeta==-1)
        beta = GetCoeffView<Scalar,typename aCoeffs::device_type>::get_view(
            -Kokkos::Details::ArithTraits<Scalar>::one(),x.dimension_1());
    }

    typedef SPMV_MV_LayoutLeft_Functor<aCoeffs,AMatrix, XVector,
                                       bCoeffs, YVector,
                                       2, 2, conjugate > OpType ;

    typename AMatrix::const_ordinal_type nrow = A.numRows();

    OpType op(alpha,A,x,beta,y,RowsPerThread<typename AMatrix::execution_space >(NNZPerRow)) ;

    // FIXME (mfh 07 Jun 2016) Shouldn't we use ordinal_type here
    // instead of int?  For example, if the number of threads is 1,
    // then this is just the number of rows.  Ditto for rows_per_team.
    // team_size is a hardware resource thing so it might legitimately
    // be int.
    const int rows_per_thread = RowsPerThread<typename AMatrix::execution_space >(NNZPerRow);
    const int team_size = Kokkos::TeamPolicy< typename AMatrix::execution_space >::team_size_recommended(op,vector_length);
    const int rows_per_team = rows_per_thread * team_size;
    const size_type nteams = (nrow+rows_per_team-1)/rows_per_team;
    Kokkos::parallel_for( Kokkos::TeamPolicy< typename AMatrix::execution_space >
       ( nteams , team_size , vector_length ) , op );

#endif // KOKKOS_FAST_COMPILE
  }
}

template<class aCoeffs,
         class AMatrix,
         class XVector,
         class bCoeffs,
         class YVector,
         int doalpha,
         int dobeta,
         bool conjugate>
static void spmv_alpha_beta_mv_transpose(aCoeffs alpha, const AMatrix& A, const XVector& x, bCoeffs beta, const YVector& y) {
  typedef typename AMatrix::ordinal_type ordinal_type;

  if (A.numRows () <= static_cast<ordinal_type> (0)) {
    return;
  }

  if (dobeta==2) {
    KokkosBlas::scal(y, beta, y);
  }
  else {
    KokkosBlas::scal(y, typename YVector::value_type (dobeta), y);
  }

  if (doalpha != 0) {
    typedef typename AMatrix::size_type size_type;

    //Impl::MV_Multiply_Check_Compatibility(betav,y,alphav,A,x,doalpha,dobeta);

    // Assuming that no row contains duplicate entries, NNZPerRow
    // cannot be more than the number of columns of the matrix.  Thus,
    // the appropriate type is ordinal_type.
    const ordinal_type NNZPerRow = static_cast<ordinal_type> (A.nnz () / A.numRows ());

    int vector_length = 1;
    while( (static_cast<ordinal_type> (vector_length*2*3) <= NNZPerRow) && (vector_length<8) ) vector_length*=2;

#ifndef KOKKOS_FAST_COMPILE // This uses templated functions on doalpha and dobeta and will produce 16 kernels

    typedef SPMV_MV_Transpose_Functor<aCoeffs,AMatrix, XVector,
                                       bCoeffs, YVector,
                                       doalpha, dobeta, conjugate > OpType ;
    OpType op(alpha,A,x,beta,y,RowsPerThread<typename AMatrix::execution_space >(NNZPerRow)) ;

    typename AMatrix::const_ordinal_type nrow = A.numRows();

    // FIXME (mfh 07 Jun 2016) Shouldn't we use ordinal_type here
    // instead of int?  For example, if the number of threads is 1,
    // then this is just the number of rows.  Ditto for rows_per_team.
    // team_size is a hardware resource thing so it might legitimately
    // be int.
    const int rows_per_thread = RowsPerThread<typename AMatrix::execution_space >(NNZPerRow);
    const int team_size = Kokkos::TeamPolicy< typename AMatrix::execution_space >::team_size_recommended(op,vector_length);
    const int rows_per_team = rows_per_thread * team_size;
    const size_type nteams = (nrow+rows_per_team-1)/rows_per_team;
    Kokkos::parallel_for( Kokkos::TeamPolicy< typename AMatrix::execution_space >
       ( nteams , team_size , vector_length ) , op );

#else // KOKKOS_FAST_COMPILE this will only instantiate one Kernel for alpha/beta

    {
      typedef typename aCoeffs::non_const_value_type Scalar;
      if(doalpha==0)
        alpha = GetCoeffView<Scalar,typename aCoeffs::device_type>::get_view(
            Kokkos::Details::ArithTraits<Scalar>::zero(),x.dimension_1());
      if(doalpha==1)
        alpha = GetCoeffView<Scalar,typename aCoeffs::device_type>::get_view(
            Kokkos::Details::ArithTraits<Scalar>::one(),x.dimension_1());
      if(doalpha==1)
        alpha = GetCoeffView<Scalar,typename aCoeffs::device_type>::get_view(
            -Kokkos::Details::ArithTraits<Scalar>::one(),x.dimension_1());
    }

    {
      typedef typename bCoeffs::non_const_value_type Scalar;
      if(dobeta==0)
        beta = GetCoeffView<Scalar,typename aCoeffs::device_type>::get_view(
            Kokkos::Details::ArithTraits<Scalar>::zero(),x.dimension_1());
      if(dobeta==1)
        beta = GetCoeffView<Scalar,typename aCoeffs::device_type>::get_view(
            Kokkos::Details::ArithTraits<Scalar>::one(),x.dimension_1());
      if(dobeta==-1)
        beta = GetCoeffView<Scalar,typename aCoeffs::device_type>::get_view(
            -Kokkos::Details::ArithTraits<Scalar>::one(),x.dimension_1());
    }

    typedef SPMV_MV_Transpose_Functor<aCoeffs,AMatrix, XVector,
                                       bCoeffs, YVector,
                                       2, 2, conjugate ,SizeType > OpType ;

    typename AMatrix::const_ordinal_type nrow = A.numRows();

    OpType op(alpha,A,x,beta,y,RowsPerThread<typename AMatrix::execution_space >(NNZPerRow)) ;

    // FIXME (mfh 07 Jun 2016) Shouldn't we use ordinal_type here
    // instead of int?  For example, if the number of threads is 1,
    // then this is just the number of rows.  Ditto for rows_per_team.
    // team_size is a hardware resource thing so it might legitimately
    // be int.
    const int rows_per_thread = RowsPerThread<typename AMatrix::execution_space >(NNZPerRow);
    const int team_size = Kokkos::TeamPolicy< typename AMatrix::execution_space >::team_size_recommended(op,vector_length);
    const int rows_per_team = rows_per_thread * team_size;
    const size_type nteams = (nrow+rows_per_team-1)/rows_per_team;
    Kokkos::parallel_for( Kokkos::TeamPolicy< typename AMatrix::execution_space >
       ( nteams , team_size , vector_length ) , op );

#endif // KOKKOS_FAST_COMPILE
  }
}

template<class aCoeffs,
         class AMatrix,
         class XVector,
         class bCoeffs,
         class YVector,
         int doalpha,
         int dobeta>
static void spmv_alpha_beta_mv(const char mode[], aCoeffs alpha, const AMatrix& A, const XVector& x, bCoeffs beta, const YVector& y) {
  if(mode[0]==NoTranspose[0]) {
    spmv_alpha_beta_mv_no_transpose<aCoeffs,AMatrix,XVector,bCoeffs,YVector,doalpha,dobeta,false>
      (alpha,A,x,beta,y);
    return;
  }
  if(mode[0]==Conjugate[0]) {
    spmv_alpha_beta_mv_no_transpose<aCoeffs,AMatrix,XVector,bCoeffs,YVector,doalpha,dobeta,true>
      (alpha,A,x,beta,y);
    return;
  }
  if(mode[0]==Transpose[0]) {
    spmv_alpha_beta_mv_transpose<aCoeffs,AMatrix,XVector,bCoeffs,YVector,doalpha,dobeta,false>
      (alpha,A,x,beta,y);
    return;
  }
  if(mode[0]==ConjugateTranspose[0]) {
    spmv_alpha_beta_mv_transpose<aCoeffs,AMatrix,XVector,bCoeffs,YVector,doalpha,dobeta,true>
      (alpha,A,x,beta,y);
    return;
  }

  Kokkos::Impl::throw_runtime_exception("Invalid Transpose Mode for KokkosSparse::spmv()");
}

template<class aCoeffs,
         class AMatrix,
         class XVector,
         class bCoeffs,
         class YVector,
         int doalpha>
void
spmv_alpha_mv (const char mode[],
               typename aCoeffs::const_value_type& alpha,
               const AMatrix& A,
               const XVector& x,
               typename bCoeffs::const_value_type& beta,
               const YVector& y)
{
  typedef typename bCoeffs::non_const_value_type Scalar;
  bCoeffs betav;

  if( beta == Kokkos::Details::ArithTraits<Scalar>::zero () ) {
    spmv_alpha_beta_mv<bCoeffs,AMatrix,XVector,bCoeffs,YVector,doalpha,0>(mode,betav,A,x,betav,y);
    return;
  }
  if( beta == Kokkos::Details::ArithTraits<Scalar>::one () ) {
    spmv_alpha_beta_mv<bCoeffs,AMatrix,XVector,bCoeffs,YVector,doalpha,1>(mode,betav,A,x,betav,y);
    return;
  }
  if( beta == -Kokkos::Details::ArithTraits<Scalar>::one () ) {
    spmv_alpha_beta_mv<bCoeffs,AMatrix,XVector,bCoeffs,YVector,doalpha,-1>(mode,betav,A,x,betav,y);
    return;
  }
  betav = GetCoeffView<Scalar,typename bCoeffs::device_type>::get_view(beta,x.dimension_1());
  spmv_alpha_beta_mv<bCoeffs,AMatrix,XVector,bCoeffs,YVector,doalpha,2>(mode,betav,A,x,betav,y);
}

#ifndef KOKKOSSPARSE_ETI_ONLY
template<class aT, class aL, class aD, class aM,
         class AT, class AO, class AD, class AM, class AS,
         class XT, class XL, class XD, class XM,
         class bT, class bL, class bD, class bM,
         class YT, class YL, class YD, class YM>
void
SPMV_MV<aT, aL, aD, aM,
        AT, AO, AD, AM, AS,
        XT, XL, XD, XM,
        bT, bL, bD, bM,
        YT, YL, YD, YM>::
spmv_mv (const char mode[],
         const aCoeffs& alpha,
         const AMatrix& A,
         const XVector& x,
         const bCoeffs& beta,
         const YVector& y)
{
  spmv_alpha_beta_mv<aCoeffs,AMatrix,XVector,bCoeffs,YVector,2,2> (mode, alpha, A, x, beta, y);
}

template<class aT, class aL, class aD, class aM,
         class AT, class AO, class AD, class AM, class AS,
         class XT, class XL, class XD, class XM,
         class bT, class bL, class bD, class bM,
         class YT, class YL, class YD, class YM>
void
SPMV_MV<aT, aL, aD, aM,
        AT, AO, AD, AM, AS,
        XT, XL, XD, XM,
        bT, bL, bD, bM,
        YT, YL, YD, YM>::
spmv_mv (const char mode[],
         const typename aCoeffs::non_const_value_type& alpha,
         const AMatrix& A,
         const XVector& x,
         const bCoeffs& beta,
         const YVector& y)
{
  typedef typename aCoeffs::non_const_value_type Scalar;

  if (alpha == Kokkos::Details::ArithTraits<Scalar>::zero ()) {
    spmv_alpha_beta_mv<bCoeffs,AMatrix,XVector,bCoeffs,YVector,0,2> (mode, beta, A, x, beta, y);
    return;
  }
  if (alpha == Kokkos::Details::ArithTraits<Scalar>::one ()) {
    spmv_alpha_beta_mv<bCoeffs,AMatrix,XVector,bCoeffs,YVector,1,2> (mode, beta, A, x, beta, y);
    return;
  }
  if (alpha == -Kokkos::Details::ArithTraits<Scalar>::one ()) {
    spmv_alpha_beta_mv<bCoeffs,AMatrix,XVector,bCoeffs,YVector,-1,2> (mode, beta, A, x, beta, y);
    return;
  }
  spmv_alpha_beta_mv<aCoeffs,AMatrix,XVector,bCoeffs,YVector,2,2> (mode,
    GetCoeffView<Scalar,typename aCoeffs::device_type>::get_view (alpha, x.dimension_1 ()),
    A, x, beta, y);
}

template<class aT, class aL, class aD, class aM,
         class AT, class AO, class AD, class AM, class AS,
         class XT, class XL, class XD, class XM,
         class bT, class bL, class bD, class bM,
         class YT, class YL, class YD, class YM>
void
SPMV_MV<aT, aL, aD, aM,
        AT, AO, AD, AM, AS,
        XT, XL, XD, XM,
        bT, bL, bD, bM,
        YT, YL, YD, YM>::
spmv_mv (const char mode[],
         const aCoeffs& alpha,
         const AMatrix& A,
         const XVector& x,
         const typename bCoeffs::non_const_value_type& beta,
         const YVector& y)
{
  typedef typename bCoeffs::non_const_value_type Scalar;

  if (beta == Kokkos::Details::ArithTraits<Scalar>::zero ()) {
    spmv_alpha_beta_mv<aCoeffs,AMatrix,XVector,aCoeffs,YVector,2,0> (mode, alpha, A, x, alpha, y);
    return;
  }
  if (beta == Kokkos::Details::ArithTraits<Scalar>::one ()) {
    spmv_alpha_beta_mv<aCoeffs,AMatrix,XVector,aCoeffs,YVector,2,1> (mode, alpha, A, x, alpha, y);
    return;
  }
  if (beta == -Kokkos::Details::ArithTraits<Scalar>::one ()) {
    spmv_alpha_beta_mv<aCoeffs,AMatrix,XVector,aCoeffs,YVector,2,-1> (mode, alpha, A, x, alpha, y);
    return;
  }
  spmv_alpha_beta_mv<aCoeffs,AMatrix,XVector,aCoeffs,YVector,2,2> (mode, alpha, A, x,
    GetCoeffView<Scalar,typename bCoeffs::device_type>::get_view (beta, x.dimension_1 ()), y);
}

template<class aT, class aL, class aD, class aM,
         class AT, class AO, class AD, class AM, class AS,
         class XT, class XL, class XD, class XM,
         class bT, class bL, class bD, class bM,
         class YT, class YL, class YD, class YM>
void
SPMV_MV<aT, aL, aD, aM,
        AT, AO, AD, AM, AS,
        XT, XL, XD, XM,
        bT, bL, bD, bM,
        YT, YL, YD, YM>::
spmv_mv (const char mode[],
         const typename aCoeffs::non_const_value_type& alpha,
         const AMatrix& A,
         const XVector& x,
         const typename bCoeffs::non_const_value_type& beta,
         const YVector& y)
{
  typedef typename aCoeffs::non_const_value_type Scalar;
  if (alpha == Kokkos::Details::ArithTraits<Scalar>::zero ()) {
    spmv_alpha_mv<aCoeffs,AMatrix,XVector,bCoeffs,YVector,0> (mode, alpha, A, x, beta, y);
    return;
  }
  if (alpha == Kokkos::Details::ArithTraits<Scalar>::one ()) {
    spmv_alpha_mv<aCoeffs,AMatrix,XVector,bCoeffs,YVector,1> (mode, alpha, A, x, beta, y);
    return;
  }
  if (alpha == -Kokkos::Details::ArithTraits<Scalar>::one ()) {
    spmv_alpha_mv<aCoeffs,AMatrix,XVector,bCoeffs,YVector,-1> (mode, alpha, A, x, beta, y);
    return;
  }
  const aCoeffs alphav = GetCoeffView<Scalar,typename aCoeffs::device_type>::get_view (alpha, x.dimension_1 ());
  spmv_mv (mode, alphav, A, x, beta, y);
}
#endif // KOKKOSSPARSE_ETI_ONLY

//
// Macro for defining (not declaring; for the declaration macro, see
// above) a full specialization of the SPMV_MV struct.  We use this
// macro in various .cpp file(s) in this directory.
//
// SCALAR_TYPE: The type of each entry in the sparse matrix
// ORDINAL_TYPE: The type of each column index in the sparse matrix
// OFFSET_TYPE: The type of each row offset in the sparse matrix
// LAYOUT_TYPE: The layout of the Kokkos::View vector arguments
//   of the sparse matrix-vector multiply
// EXEC_SPACE_TYPE: The execution space type
// MEM_SPACE_TYPE: The memory space type
//
#define KOKKOSSPARSE_IMPL_SPMV_MV_DEF( SCALAR_TYPE, ORDINAL_TYPE, OFFSET_TYPE, LAYOUT_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE ) \
void \
SPMV_MV<const SCALAR_TYPE*, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
        const SCALAR_TYPE, \
        ORDINAL_TYPE, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
        OFFSET_TYPE, \
        const SCALAR_TYPE**, \
        LAYOUT_TYPE, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess>, \
        const SCALAR_TYPE*, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
        SCALAR_TYPE**, \
        LAYOUT_TYPE, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> >:: \
spmv_mv (const char mode[], \
         const aCoeffs& alpha, \
         const AMatrix& A, \
         const XVector& x, \
         const bCoeffs& beta, \
         const YVector& y) \
{ \
  spmv_alpha_beta_mv<aCoeffs,AMatrix,XVector,bCoeffs,YVector,2,2> (mode, alpha, A, x, beta, y); \
} \
  \
void \
SPMV_MV<const SCALAR_TYPE*, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> , \
        const SCALAR_TYPE, \
        ORDINAL_TYPE, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
        OFFSET_TYPE, \
        const SCALAR_TYPE**, \
        LAYOUT_TYPE, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> , \
        const SCALAR_TYPE*, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> , \
        SCALAR_TYPE**, \
        LAYOUT_TYPE, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> >:: \
spmv_mv (const char mode[], \
         const aCoeffs::non_const_value_type& alpha, \
         const AMatrix& A, \
         const XVector& x, \
         const bCoeffs& beta, \
         const YVector& y) \
{ \
  typedef aCoeffs::non_const_value_type Scalar; \
  if (alpha == Kokkos::Details::ArithTraits<Scalar>::zero ()) { \
    spmv_alpha_beta_mv<bCoeffs,AMatrix,XVector,bCoeffs,YVector,0,2> (mode, beta, A, x, beta, y); \
    return; \
  } \
  if (alpha == Kokkos::Details::ArithTraits<Scalar>::one ()) { \
    spmv_alpha_beta_mv<bCoeffs,AMatrix,XVector,bCoeffs,YVector,1,2> (mode, beta, A, x, beta, y); \
    return; \
  } \
  if (alpha == -Kokkos::Details::ArithTraits<Scalar>::one ()) { \
    spmv_alpha_beta_mv<bCoeffs,AMatrix,XVector,bCoeffs,YVector,-1,2> (mode, beta, A, x, beta, y); \
    return; \
  } \
  spmv_alpha_beta_mv<aCoeffs,AMatrix,XVector,bCoeffs,YVector,2,2> (mode, \
    GetCoeffView<Scalar,aCoeffs::device_type>::get_view (alpha, x.dimension_1 ()), \
    A, x, beta, y); \
} \
  \
void \
SPMV_MV<const SCALAR_TYPE*, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> , \
        const SCALAR_TYPE, \
        ORDINAL_TYPE, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
        OFFSET_TYPE, \
        const SCALAR_TYPE**, \
        LAYOUT_TYPE, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> , \
        const SCALAR_TYPE*, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> , \
        SCALAR_TYPE**, \
        LAYOUT_TYPE, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> >:: \
spmv_mv (const char mode[], \
         const aCoeffs& alpha, \
         const AMatrix& A, \
         const XVector& x, \
         const bCoeffs::non_const_value_type& beta, \
         const YVector& y) \
{ \
  typedef bCoeffs::non_const_value_type Scalar; \
  if (beta == Kokkos::Details::ArithTraits<Scalar>::zero ()) { \
    spmv_alpha_beta_mv<aCoeffs,AMatrix,XVector,aCoeffs,YVector,2,0> (mode, alpha, A, x, alpha, y); \
    return; \
  } \
  if (beta == Kokkos::Details::ArithTraits<Scalar>::one ()) { \
    spmv_alpha_beta_mv<aCoeffs,AMatrix,XVector,aCoeffs,YVector,2,1> (mode, alpha, A, x, alpha, y); \
    return; \
  } \
  if (beta == -Kokkos::Details::ArithTraits<Scalar>::one ()) { \
    spmv_alpha_beta_mv<aCoeffs,AMatrix,XVector,aCoeffs,YVector,2,-1> (mode, alpha, A, x, alpha, y); \
    return; \
  } \
  spmv_alpha_beta_mv<aCoeffs,AMatrix,XVector,aCoeffs,YVector,2,2> (mode, alpha, A, x, \
    GetCoeffView<Scalar,bCoeffs::device_type>::get_view (beta, x.dimension_1 ()), y); \
} \
  \
void \
SPMV_MV<const SCALAR_TYPE*, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> , \
        const SCALAR_TYPE, \
        ORDINAL_TYPE, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
        OFFSET_TYPE, \
        const SCALAR_TYPE**, \
        LAYOUT_TYPE, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> , \
        const SCALAR_TYPE*, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> , \
        SCALAR_TYPE**, \
        LAYOUT_TYPE, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> >:: \
spmv_mv (const char mode[], \
         const aCoeffs::non_const_value_type& alpha, \
         const AMatrix& A, \
         const XVector& x, \
         const bCoeffs::non_const_value_type& beta, \
         const YVector& y) \
{ \
  typedef aCoeffs::non_const_value_type Scalar; \
  if (alpha == Kokkos::Details::ArithTraits<Scalar>::zero ()) { \
    spmv_alpha_mv<aCoeffs,AMatrix,XVector,bCoeffs,YVector,0> (mode, alpha, A, x, beta, y); \
    return; \
  } \
  if (alpha == Kokkos::Details::ArithTraits<Scalar>::one ()) { \
    spmv_alpha_mv<aCoeffs,AMatrix,XVector,bCoeffs,YVector,1> (mode, alpha, A, x, beta, y); \
    return; \
  } \
  if (alpha == -Kokkos::Details::ArithTraits<Scalar>::one ()) { \
    spmv_alpha_mv<aCoeffs,AMatrix,XVector,bCoeffs,YVector,-1> (mode, alpha, A, x, beta, y); \
    return; \
  } \
  const aCoeffs alphav = GetCoeffView<Scalar,aCoeffs::device_type>::get_view (alpha, x.dimension_1 ()); \
  spmv_mv (mode, alphav, A, x, beta, y); \
}


//
// Macro for defining (not declaring; for the declaration macro, see
// above) a full specialization of the SPMV_MV struct, with the
// default OFFSET_TYPE and LAYOUT_TYPE = Kokkos::LayoutLeft (Tpetra's
// layout).  We use this macro in various .cpp file(s) in this
// directory.
//
// SCALAR_TYPE: The type of each entry in the sparse matrix
// ORDINAL_TYPE: The type of each column index in the sparse matrix
// EXEC_SPACE_TYPE: The execution space type
// MEM_SPACE_TYPE: The memory space type
//
#define KOKKOSSPARSE_IMPL_SPMV_MV_DEFAULTS_DEF( SCALAR_TYPE, ORDINAL_TYPE, EXEC_SPACE_TYPE, MEM_SPACE_TYPE ) \
void \
SPMV_MV<const SCALAR_TYPE*, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
        const SCALAR_TYPE, \
        ORDINAL_TYPE, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
        Kokkos::StaticCrsGraph<ORDINAL_TYPE, Kokkos::LayoutLeft, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE> >::size_type, \
        const SCALAR_TYPE**, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess>, \
        const SCALAR_TYPE*, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
        SCALAR_TYPE**, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> >:: \
spmv_mv (const char mode[], \
         const aCoeffs& alpha, \
         const AMatrix& A, \
         const XVector& x, \
         const bCoeffs& beta, \
         const YVector& y) \
{ \
  spmv_alpha_beta_mv<aCoeffs,AMatrix,XVector,bCoeffs,YVector,2,2> (mode, alpha, A, x, beta, y); \
} \
  \
void \
SPMV_MV<const SCALAR_TYPE*, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> , \
        const SCALAR_TYPE, \
        ORDINAL_TYPE, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
        Kokkos::StaticCrsGraph<ORDINAL_TYPE, Kokkos::LayoutLeft, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE> >::size_type, \
        const SCALAR_TYPE**, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> , \
        const SCALAR_TYPE*, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> , \
        SCALAR_TYPE**, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> >:: \
spmv_mv (const char mode[], \
         const aCoeffs::non_const_value_type& alpha, \
         const AMatrix& A, \
         const XVector& x, \
         const bCoeffs& beta, \
         const YVector& y) \
{ \
  typedef aCoeffs::non_const_value_type Scalar; \
  if (alpha == Kokkos::Details::ArithTraits<Scalar>::zero ()) { \
    spmv_alpha_beta_mv<bCoeffs,AMatrix,XVector,bCoeffs,YVector,0,2> (mode, beta, A, x, beta, y); \
    return; \
  } \
  if (alpha == Kokkos::Details::ArithTraits<Scalar>::one ()) { \
    spmv_alpha_beta_mv<bCoeffs,AMatrix,XVector,bCoeffs,YVector,1,2> (mode, beta, A, x, beta, y); \
    return; \
  } \
  if (alpha == -Kokkos::Details::ArithTraits<Scalar>::one ()) { \
    spmv_alpha_beta_mv<bCoeffs,AMatrix,XVector,bCoeffs,YVector,-1,2> (mode, beta, A, x, beta, y); \
    return; \
  } \
  spmv_alpha_beta_mv<aCoeffs,AMatrix,XVector,bCoeffs,YVector,2,2> (mode, \
    GetCoeffView<Scalar,aCoeffs::device_type>::get_view (alpha, x.dimension_1 ()), \
    A, x, beta, y); \
} \
  \
void \
SPMV_MV<const SCALAR_TYPE*, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> , \
        const SCALAR_TYPE, \
        ORDINAL_TYPE, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
        Kokkos::StaticCrsGraph<ORDINAL_TYPE, Kokkos::LayoutLeft, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE> >::size_type, \
        const SCALAR_TYPE**, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> , \
        const SCALAR_TYPE*, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> , \
        SCALAR_TYPE**, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> >:: \
spmv_mv (const char mode[], \
         const aCoeffs& alpha, \
         const AMatrix& A, \
         const XVector& x, \
         const bCoeffs::non_const_value_type& beta, \
         const YVector& y) \
{ \
  typedef bCoeffs::non_const_value_type Scalar; \
  if (beta == Kokkos::Details::ArithTraits<Scalar>::zero ()) { \
    spmv_alpha_beta_mv<aCoeffs,AMatrix,XVector,aCoeffs,YVector,2,0> (mode, alpha, A, x, alpha, y); \
    return; \
  } \
  if (beta == Kokkos::Details::ArithTraits<Scalar>::one ()) { \
    spmv_alpha_beta_mv<aCoeffs,AMatrix,XVector,aCoeffs,YVector,2,1> (mode, alpha, A, x, alpha, y); \
    return; \
  } \
  if (beta == -Kokkos::Details::ArithTraits<Scalar>::one ()) { \
    spmv_alpha_beta_mv<aCoeffs,AMatrix,XVector,aCoeffs,YVector,2,-1> (mode, alpha, A, x, alpha, y); \
    return; \
  } \
  spmv_alpha_beta_mv<aCoeffs,AMatrix,XVector,aCoeffs,YVector,2,2> (mode, alpha, A, x, \
    GetCoeffView<Scalar,bCoeffs::device_type>::get_view (beta, x.dimension_1 ()), y); \
} \
  \
void \
SPMV_MV<const SCALAR_TYPE*, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> , \
        const SCALAR_TYPE, \
        ORDINAL_TYPE, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged>, \
        Kokkos::StaticCrsGraph<ORDINAL_TYPE, Kokkos::LayoutLeft, Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE> >::size_type, \
        const SCALAR_TYPE**, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged|Kokkos::RandomAccess> , \
        const SCALAR_TYPE*, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> , \
        SCALAR_TYPE**, \
        Kokkos::LayoutLeft, \
        Kokkos::Device<EXEC_SPACE_TYPE, MEM_SPACE_TYPE>, \
        Kokkos::MemoryTraits<Kokkos::Unmanaged> >:: \
spmv_mv (const char mode[], \
         const aCoeffs::non_const_value_type& alpha, \
         const AMatrix& A, \
         const XVector& x, \
         const bCoeffs::non_const_value_type& beta, \
         const YVector& y) \
{ \
  typedef aCoeffs::non_const_value_type Scalar; \
  if (alpha == Kokkos::Details::ArithTraits<Scalar>::zero ()) { \
    spmv_alpha_mv<aCoeffs,AMatrix,XVector,bCoeffs,YVector,0> (mode, alpha, A, x, beta, y); \
    return; \
  } \
  if (alpha == Kokkos::Details::ArithTraits<Scalar>::one ()) { \
    spmv_alpha_mv<aCoeffs,AMatrix,XVector,bCoeffs,YVector,1> (mode, alpha, A, x, beta, y); \
    return; \
  } \
  if (alpha == -Kokkos::Details::ArithTraits<Scalar>::one ()) { \
    spmv_alpha_mv<aCoeffs,AMatrix,XVector,bCoeffs,YVector,-1> (mode, alpha, A, x, beta, y); \
    return; \
  } \
  const aCoeffs alphav = GetCoeffView<Scalar,aCoeffs::device_type>::get_view (alpha, x.dimension_1 ()); \
  spmv_mv (mode, alphav, A, x, beta, y); \
}

} // namespace Impl
} // namespace KokkosSparse

#endif // KOKKOS_SPARSE_IMPL_SPMV_DEF_HPP_
