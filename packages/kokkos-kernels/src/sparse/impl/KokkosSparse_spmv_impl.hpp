/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOSSPARSE_IMPL_SPMV_DEF_HPP_
#define KOKKOSSPARSE_IMPL_SPMV_DEF_HPP_

#include "KokkosKernels_Controls.hpp"
#include "Kokkos_InnerProductSpaceTraits.hpp"
#include "KokkosBlas1_scal.hpp"
#include "KokkosSparse_CrsMatrix.hpp"
#include "KokkosSparse_spmv_impl_omp.hpp"

namespace KokkosSparse {
namespace Impl {

template<class InputType, class DeviceType>
struct GetCoeffView {
  typedef Kokkos::View<InputType*,Kokkos::LayoutLeft,DeviceType> view_type;
  typedef Kokkos::View<typename view_type::non_const_value_type*,
                       Kokkos::LayoutLeft,DeviceType> non_const_view_type;
  static non_const_view_type get_view(const InputType in, const int size) {
    non_const_view_type aview("CoeffView",size);
    if(size>0)
      Kokkos::deep_copy(aview,in);
    return aview;
  }
};

template<class IT, class IL, class ID, class IM, class IS, class DeviceType>
struct GetCoeffView<Kokkos::View<IT*,IL,ID,IM,IS>,DeviceType> {
  typedef Kokkos::View<IT*,IL,ID,IM,IS> view_type;
  static Kokkos::View<IT*,IL,ID,IM,IS> get_view(const Kokkos::View<IT*,IL,ID,IM,IS>& in, int size) {
    return in;
  }
};


// This TransposeFunctor is functional, but not necessarily performant.
template<class AMatrix,
         class XVector,
         class YVector,
         int dobeta,
         bool conjugate>
struct SPMV_Transpose_Functor {
  typedef typename AMatrix::execution_space            execution_space;
  typedef typename AMatrix::non_const_ordinal_type     ordinal_type;
  typedef typename AMatrix::non_const_value_type       value_type;
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type            team_member;
  typedef Kokkos::Details::ArithTraits<value_type>     ATV;
  typedef typename YVector::non_const_value_type       coefficient_type;
  typedef typename YVector::non_const_value_type       y_value_type;

  const coefficient_type alpha;
  AMatrix m_A;
  XVector m_x;
  const coefficient_type beta;
  YVector m_y;
  const ordinal_type rows_per_thread;

  SPMV_Transpose_Functor (const coefficient_type& alpha_,
                          const AMatrix& m_A_,
                          const XVector& m_x_,
                          const coefficient_type& beta_,
                          const YVector& m_y_,
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
           iEntry += static_cast<ordinal_type> (blockDim.x))
#else
      for (ordinal_type iEntry = 0;
           iEntry < row_length;
           iEntry ++)
#endif
      {
        const value_type val = conjugate ?
          ATV::conj (row.value(iEntry)) :
          row.value(iEntry);
        const ordinal_type ind = row.colidx(iEntry);

        Kokkos::atomic_add (&m_y(ind), static_cast<y_value_type> (alpha * val * m_x(iRow)));
      }
    }
  }
};

template<class AMatrix,
         class XVector,
         class YVector,
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

  const ordinal_type rows_per_team;

  SPMV_Functor (const value_type alpha_,
                const AMatrix m_A_,
                const XVector m_x_,
                const value_type beta_,
                const YVector m_y_,
                const int rows_per_team_) :
     alpha (alpha_), m_A (m_A_), m_x (m_x_),
     beta (beta_), m_y (m_y_),
     rows_per_team (rows_per_team_)
  {
    static_assert (static_cast<int> (XVector::rank) == 1,
                   "XVector must be a rank 1 View.");
    static_assert (static_cast<int> (YVector::rank) == 1,
                   "YVector must be a rank 1 View.");
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const team_member& dev) const
  {
    typedef typename YVector::non_const_value_type y_value_type;

    Kokkos::parallel_for(Kokkos::TeamThreadRange(dev,0,rows_per_team), [&] (const ordinal_type& loop) {

      const ordinal_type iRow = static_cast<ordinal_type> ( dev.league_rank() ) * rows_per_team + loop;
      if (iRow >= m_A.numRows ()) {
        return;
      }
      const KokkosSparse::SparseRowViewConst<AMatrix> row = m_A.rowConst(iRow);
      const ordinal_type row_length = static_cast<ordinal_type> (row.length);
      y_value_type sum = 0;

      Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(dev,row_length), [&] (const ordinal_type& iEntry, y_value_type& lsum) {
        const value_type val = conjugate ?
                ATV::conj (row.value(iEntry)) :
                row.value(iEntry);
        lsum += val * m_x(row.colidx(iEntry));
      },sum);

      Kokkos::single(Kokkos::PerThread(dev), [&] () {
        sum *= alpha;

        if (dobeta == 0) {
          m_y(iRow) = sum ;
        } else {
          m_y(iRow) = beta * m_y(iRow) + sum;
        }
      });
    });
  }
};

template<class execution_space>
int64_t spmv_launch_parameters(int64_t numRows, int64_t nnz, int64_t rows_per_thread, int& team_size, int& vector_length) {
  int64_t rows_per_team;
  int64_t nnz_per_row = nnz/numRows;

  if(nnz_per_row < 1) nnz_per_row = 1;

  if(vector_length < 1) {
    vector_length = 1;
    while(vector_length<32 && vector_length*6 < nnz_per_row)
      vector_length*=2;
  }

  // Determine rows per thread
  if(rows_per_thread < 1) {
    #ifdef KOKKOS_ENABLE_CUDA
    if(std::is_same<Kokkos::Cuda,execution_space>::value)
      rows_per_thread = 1;
    else
    #endif
    {
      if(nnz_per_row < 20 && nnz > 5000000 ) {
        rows_per_thread = 256;
      } else
        rows_per_thread = 64;
    }
  }

  #ifdef KOKKOS_ENABLE_CUDA
  if(team_size < 1) {
    if(std::is_same<Kokkos::Cuda,execution_space>::value)
    { team_size = 256/vector_length; }
    else
    { team_size = 1; }
  }
  #endif

  rows_per_team = rows_per_thread * team_size;

  if(rows_per_team < 0) {
    int64_t nnz_per_team = 4096;
    int64_t conc = execution_space::concurrency();
    while((conc * nnz_per_team * 4> nnz)&&(nnz_per_team>256)) nnz_per_team/=2;
    rows_per_team = (nnz_per_team+nnz_per_row - 1)/nnz_per_row;
  }


  return rows_per_team;
}

template<class AMatrix,
         class XVector,
         class YVector,
         int dobeta,
         bool conjugate>
static void
spmv_beta_no_transpose (const KokkosKernels::Experimental::Controls& controls,
			typename YVector::const_value_type& alpha,
			const AMatrix& A,
			const XVector& x,
			typename YVector::const_value_type& beta,
			const YVector& y)
{
  typedef typename AMatrix::ordinal_type ordinal_type;
  typedef typename AMatrix::execution_space execution_space;

  if (A.numRows () <= static_cast<ordinal_type> (0)) {
    return;
  }

  #ifdef KOKKOS_ENABLE_OPENMP
  if((std::is_same<execution_space,Kokkos::OpenMP>::value) &&
     (std::is_same<typename std::remove_cv<typename AMatrix::value_type>::type,double>::value) &&
     (std::is_same<typename XVector::non_const_value_type,double>::value) &&
     (std::is_same<typename YVector::non_const_value_type,double>::value) &&
     ((int) A.graph.row_block_offsets.extent(0) == (int) omp_get_max_threads()+1) &&
     (((uintptr_t)(const void*)(x.data())%64)==0) && (((uintptr_t)(const void*)(y.data())%64)==0)
     ) {
    spmv_raw_openmp_no_transpose<AMatrix,XVector,YVector>(alpha,A,x,beta,y);
    return;
  }
  #endif
  int team_size = -1;
  int vector_length = -1;
  int64_t rows_per_thread = -1;

  // Note on 03/24/20, lbv: We can use the controls
  // here to allow the user to pass in some tunning
  // parameters.
  if(controls.isParameter("team size"))       {team_size       = std::stoi(controls.getParameter("team size"));}
  if(controls.isParameter("vector length"))   {vector_length   = std::stoi(controls.getParameter("vector length"));}
  if(controls.isParameter("rows per thread")) {rows_per_thread = std::stoll(controls.getParameter("rows per thread"));}

  bool use_dynamic_schedule = false; // Forces the use of a dynamic schedule
  bool use_static_schedule  = false; // Forces the use of a static schedule
  if(controls.isParameter("schedule")) {
    if(controls.getParameter("schedule") == "dynamic") {
      use_dynamic_schedule = true;
    } else if(controls.getParameter("schedule") == "static") {
      use_static_schedule  = true;
    }
  }

  int64_t rows_per_team = spmv_launch_parameters<execution_space>(A.numRows(),A.nnz(),rows_per_thread,team_size,vector_length);
  int64_t worksets = (y.extent(0)+rows_per_team-1)/rows_per_team;

  SPMV_Functor<AMatrix,XVector,YVector,dobeta,conjugate> func (alpha,A,x,beta,y,rows_per_team);

  if(((A.nnz()>10000000) || use_dynamic_schedule) && !use_static_schedule) {
    Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic> > policy(1,1);
    if(team_size<0)
      policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic> >(worksets,Kokkos::AUTO,vector_length);
    else
      policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Dynamic> >(worksets,team_size,vector_length);
    Kokkos::parallel_for("KokkosSparse::spmv<NoTranspose,Dynamic>",policy,func);
  } else {
    Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static> > policy(1,1);
    if(team_size<0)
      policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static> >(worksets,Kokkos::AUTO,vector_length);
    else
      policy = Kokkos::TeamPolicy<execution_space, Kokkos::Schedule<Kokkos::Static> >(worksets,team_size,vector_length);
    Kokkos::parallel_for("KokkosSparse::spmv<NoTranspose,Static>",policy,func);
  }
}

template<class AMatrix,
         class XVector,
         class YVector,
         int dobeta,
         bool conjugate>
static void
spmv_beta_transpose (typename YVector::const_value_type& alpha,
                           const AMatrix& A,
                           const XVector& x,
                           typename YVector::const_value_type& beta,
                           const YVector& y)
{
  typedef typename AMatrix::ordinal_type ordinal_type;

  if (A.numRows () <= static_cast<ordinal_type> (0)) {
    return;
  }

  // We need to scale y first ("scaling" by zero just means filling
  // with zeros), since the functor works by atomic-adding into y.
  if (dobeta != 1) {
    KokkosBlas::scal (y, beta, y);
  }

  typedef typename AMatrix::size_type size_type;

  // Assuming that no row contains duplicate entries, NNZPerRow
  // cannot be more than the number of columns of the matrix.  Thus,
  // the appropriate type is ordinal_type.
  const ordinal_type NNZPerRow = static_cast<ordinal_type> (A.nnz () / A.numRows ());

  int vector_length = 1;
  while( (static_cast<ordinal_type> (vector_length*2*3) <= NNZPerRow) && (vector_length<32) ) vector_length*=2;

  typedef SPMV_Transpose_Functor<AMatrix, XVector, YVector, dobeta, conjugate> OpType;

  typename AMatrix::const_ordinal_type nrow = A.numRows();

  OpType op (alpha, A, x, beta, y, RowsPerThread<typename AMatrix::execution_space> (NNZPerRow));

  const int rows_per_thread = RowsPerThread<typename AMatrix::execution_space > (NNZPerRow);
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
  const int team_size = Kokkos::TeamPolicy<typename AMatrix::execution_space>::team_size_recommended (op, vector_length);
#else
  const int team_size = Kokkos::TeamPolicy<typename AMatrix::execution_space>(rows_per_thread, Kokkos::AUTO, vector_length).team_size_recommended(op, Kokkos::ParallelForTag());
#endif
  const int rows_per_team = rows_per_thread * team_size;
  const size_type nteams = (nrow+rows_per_team-1)/rows_per_team;
  Kokkos::parallel_for("KokkosSparse::spmv<Transpose>", Kokkos::TeamPolicy< typename AMatrix::execution_space >
     ( nteams , team_size , vector_length ) , op );

}

template<class AMatrix,
         class XVector,
         class YVector,
         int dobeta>
static void
spmv_beta (const KokkosKernels::Experimental::Controls& controls,
	   const char mode[],
	   typename YVector::const_value_type& alpha,
	   const AMatrix& A,
	   const XVector& x,
	   typename YVector::const_value_type& beta,
	   const YVector& y)
{
  if (mode[0] == NoTranspose[0]) {
    spmv_beta_no_transpose<AMatrix,XVector,YVector,dobeta,false>
      (controls,alpha,A,x,beta,y);
  }
  else if (mode[0] == Conjugate[0]) {
    spmv_beta_no_transpose<AMatrix,XVector,YVector,dobeta,true>
      (controls,alpha,A,x,beta,y);
  }
  else if (mode[0]==Transpose[0]) {
    spmv_beta_transpose<AMatrix,XVector,YVector,dobeta,false>
      (alpha,A,x,beta,y);
  }
  else if(mode[0]==ConjugateTranspose[0]) {
    spmv_beta_transpose<AMatrix,XVector,YVector,dobeta,true>
      (alpha,A,x,beta,y);
  }
  else {
    Kokkos::Impl::throw_runtime_exception("Invalid Transpose Mode for KokkosSparse::spmv()");
  }
}


// Functor for implementing transpose and conjugate transpose sparse
// matrix-vector multiply with multivector (2-D View) input and
// output.  This functor works, but is not necessarily performant.
template<class AMatrix,
         class XVector,
         class YVector,
         int doalpha,
         int dobeta,
         bool conjugate>
struct SPMV_MV_Transpose_Functor {
  typedef typename AMatrix::execution_space            execution_space;
  typedef typename AMatrix::non_const_ordinal_type     ordinal_type;
  typedef typename AMatrix::non_const_value_type       A_value_type;
  typedef typename YVector::non_const_value_type       y_value_type;
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type            team_member;
  typedef typename YVector::non_const_value_type       coefficient_type;

  const coefficient_type alpha;
  AMatrix m_A;
  XVector m_x;
  const coefficient_type beta;
  YVector m_y;

  const ordinal_type n;
  const ordinal_type rows_per_thread;

  SPMV_MV_Transpose_Functor (const coefficient_type& alpha_,
                             const AMatrix& m_A_,
                             const XVector& m_x_,
                             const coefficient_type& beta_,
                             const YVector& m_y_,
                             const ordinal_type rows_per_thread_) :
    alpha (alpha_),
    m_A (m_A_), m_x (m_x_), beta (beta_), m_y (m_y_), n (m_x_.extent(1)),
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
           iEntry += static_cast<ordinal_type> (blockDim.x))
#else
      for (ordinal_type iEntry = 0;
           iEntry < row_length;
           iEntry ++)
#endif
      {
        const A_value_type val = conjugate ?
          Kokkos::Details::ArithTraits<A_value_type>::conj (row.value(iEntry)) :
          row.value(iEntry);
        const ordinal_type ind = row.colidx(iEntry);

        if (doalpha != 1) {
          #ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
          #pragma unroll
          #endif
          for (ordinal_type k = 0; k < n; ++k) {
            Kokkos::atomic_add (&m_y(ind,k),
                                static_cast<y_value_type> (alpha * val * m_x(iRow, k)));
          }
        } else {
          #ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
          #pragma unroll
          #endif
          for (ordinal_type k = 0; k < n; ++k) {
            Kokkos::atomic_add (&m_y(ind,k),
                                static_cast<y_value_type> (val * m_x(iRow, k)));
          }
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
struct SPMV_MV_LayoutLeft_Functor {
  typedef typename AMatrix::execution_space            execution_space;
  typedef typename AMatrix::non_const_ordinal_type               ordinal_type;
  typedef typename AMatrix::non_const_value_type       A_value_type;
  typedef typename YVector::non_const_value_type       y_value_type;
  typedef typename Kokkos::TeamPolicy<execution_space> team_policy;
  typedef typename team_policy::member_type            team_member;
  typedef typename YVector::non_const_value_type       coefficient_type;

  const coefficient_type alpha;
  AMatrix m_A;
  XVector m_x;
  const coefficient_type beta;
  YVector m_y;
  //! The number of columns in the input and output MultiVectors.
  ordinal_type n;
  ordinal_type rows_per_thread;

  SPMV_MV_LayoutLeft_Functor (const coefficient_type& alpha_,
                              const AMatrix& m_A_,
                              const XVector& m_x_,
                              const coefficient_type& beta_,
                              const YVector& m_y_,
                              const ordinal_type rows_per_thread_) :
    alpha (alpha_),
    m_A (m_A_), m_x (m_x_), beta (beta_), m_y (m_y_), n (m_x_.extent(1)),
    rows_per_thread (rows_per_thread_)
  {}

  template<int UNROLL>
  KOKKOS_INLINE_FUNCTION void
  strip_mine (const team_member& /* dev */, const ordinal_type& iRow, const ordinal_type& kk) const
  {
    y_value_type sum[UNROLL];

#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
    for (int k = 0; k < UNROLL; ++k) {
      sum[k] = Kokkos::Details::ArithTraits<y_value_type>::zero ();
    }

    const auto row = m_A.rowConst (iRow);

    // The correct type of iEntry is ordinal_type, the type of the
    // number of columns in the (local) matrix.  This is because we
    // assume either that rows have no duplicate entries, or that rows
    // never have enough duplicate entries to overflow ordinal_type.

#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_LOOPCOUNT
#pragma loop count (15)
#endif
#ifdef __CUDA_ARCH__
        for (ordinal_type iEntry = static_cast<ordinal_type> (threadIdx.x);
             iEntry < row.length;
             iEntry += static_cast<ordinal_type> (blockDim.x))
#else
        for (ordinal_type iEntry = 0;
             iEntry < row.length;
             iEntry ++)
#endif
      {
      const A_value_type val = conjugate ?
        Kokkos::Details::ArithTraits<A_value_type>::conj (row.value(iEntry)) :
        row.value(iEntry);
      const ordinal_type ind = row.colidx(iEntry);

#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
      for (int k = 0; k < UNROLL; ++k) {
        sum[k] += val * m_x(ind, kk + k);
      }
    }

    if (doalpha == -1) {
      for (int ii=0; ii < UNROLL; ++ii) {
        y_value_type sumt = sum[ii];
#if defined(__CUDA_ARCH__) && defined(KOKKOS_ENABLE_CUDA)
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
#endif // defined(__CUDA_ARCH__) && defined(KOKKOS_ENABLE_CUDA)
        sum[ii] = -sumt;
      }
    }
    else {
      for (int ii=0; ii < UNROLL; ++ii) {
        y_value_type sumt = sum[ii];
#if defined(__CUDA_ARCH__) && defined(KOKKOS_ENABLE_CUDA)
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
#endif // defined(__CUDA_ARCH__) && defined(KOKKOS_ENABLE_CUDA)
        sum[ii] = sumt;
      }
    }

#if defined(__CUDA_ARCH__) && defined(KOKKOS_ENABLE_CUDA)
    if (threadIdx.x==0)
#else
    if (true)
#endif // defined(__CUDA_ARCH__) && defined(KOKKOS_ENABLE_CUDA)
    {
      if (doalpha * doalpha != 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          sum[k] *= alpha;
        }
      }

      if (dobeta == 0) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(iRow, kk + k) = sum[k];
        }
      } else if (dobeta == 1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(iRow, kk + k) += sum[k];
        }
      } else if (dobeta == -1) {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(iRow, kk + k) = -m_y(iRow, kk + k) +  sum[k];
        }
      } else {
#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
        for (int k = 0; k < UNROLL; ++k) {
          m_y(iRow, kk + k) = beta * m_y(iRow, kk + k) + sum[k];
        }
      }
    }
  }

  KOKKOS_INLINE_FUNCTION void
  strip_mine_1 (const team_member& /* dev */, const ordinal_type& iRow) const
  {
    y_value_type sum = Kokkos::Details::ArithTraits<y_value_type>::zero ();

    const auto row = m_A.rowConst (iRow);

    // The correct type of iEntry is ordinal_type, the type of the
    // number of columns in the (local) matrix.  This is because we
    // assume either that rows have no duplicate entries, or that rows
    // never have enough duplicate entries to overflow ordinal_type.

#ifdef KOKKOS_ENABLE_PRAGMA_IVDEP
#pragma ivdep
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_UNROLL
#pragma unroll
#endif
#ifdef KOKKOS_ENABLE_PRAGMA_LOOPCOUNT
#pragma loop count (15)
#endif
#ifdef __CUDA_ARCH__
    for (ordinal_type iEntry = static_cast<ordinal_type> (threadIdx.x);
         iEntry < row.length;
         iEntry += static_cast<ordinal_type> (blockDim.x))
#else
    for (ordinal_type iEntry = 0;
         iEntry < row.length;
         iEntry ++)
#endif
    {
      const A_value_type val = conjugate ?
          Kokkos::Details::ArithTraits<A_value_type>::conj (row.value(iEntry)) :
          row.value(iEntry);
      sum += val * m_x(row.colidx(iEntry),0);
    }
#if defined(__CUDA_ARCH__) && defined(KOKKOS_ENABLE_CUDA)
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
#endif // defined(__CUDA_ARCH__) && defined(KOKKOS_ENABLE_CUDA)

#if defined(__CUDA_ARCH__) && defined(KOKKOS_ENABLE_CUDA)
    if (threadIdx.x==0)
#else
    if (true)
#endif // defined(__CUDA_ARCH__) && defined(KOKKOS_ENABLE_CUDA)
    {
      if (doalpha == -1) {
        sum = -sum;
      } else if (doalpha * doalpha != 1) {
        sum *= alpha;
      }

      if (dobeta == 0) {
        m_y(iRow, 0) = sum ;
      } else if (dobeta == 1) {
        m_y(iRow, 0) += sum ;
      } else if (dobeta == -1) {
        m_y(iRow, 0) = -m_y(iRow, 0) +  sum;
      } else {
        m_y(iRow, 0) = beta * m_y(iRow, 0) + sum;
      }
    }
  }


  KOKKOS_INLINE_FUNCTION void
  operator() (const team_member& dev) const
  {
    for (ordinal_type loop = 0; loop < rows_per_thread; ++loop) {

      // iRow indexes over (local) rows of the matrix, so its correct
      // type is ordinal_type.

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


template<class AMatrix,
         class XVector,
         class YVector,
         int doalpha,
         int dobeta,
         bool conjugate>
static void
spmv_alpha_beta_mv_no_transpose (const typename YVector::non_const_value_type& alpha,
                                 const AMatrix& A,
                                 const XVector& x,
                                 const typename YVector::non_const_value_type& beta,
                                 const YVector& y)
{
  typedef typename AMatrix::ordinal_type ordinal_type;

  if (A.numRows () <= static_cast<ordinal_type> (0)) {
    return;
  }
  if (doalpha == 0) {
    if (dobeta != 1) {
      KokkosBlas::scal (y, beta, y);
    }
    return;
  }
  else {
    typedef typename AMatrix::size_type size_type;

    // Assuming that no row contains duplicate entries, NNZPerRow
    // cannot be more than the number of columns of the matrix.  Thus,
    // the appropriate type is ordinal_type.
    const ordinal_type NNZPerRow = static_cast<ordinal_type> (A.nnz () / A.numRows ());

    int vector_length = 1;
    while( (static_cast<ordinal_type> (vector_length*2*3) <= NNZPerRow) && (vector_length<8) ) vector_length*=2;

#ifndef KOKKOS_FAST_COMPILE // This uses templated functions on doalpha and dobeta and will produce 16 kernels

    typedef SPMV_MV_LayoutLeft_Functor<AMatrix, XVector, YVector,
                                       doalpha, dobeta, conjugate> OpType;
    OpType op (alpha, A, x, beta, y, RowsPerThread<typename AMatrix::execution_space> (NNZPerRow));

    typename AMatrix::const_ordinal_type nrow = A.numRows();

    // FIXME (mfh 07 Jun 2016) Shouldn't we use ordinal_type here
    // instead of int?  For example, if the number of threads is 1,
    // then this is just the number of rows.  Ditto for rows_per_team.
    // team_size is a hardware resource thing so it might legitimately
    // be int.
    const int rows_per_thread = RowsPerThread<typename AMatrix::execution_space >(NNZPerRow);
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    const int team_size = Kokkos::TeamPolicy< typename AMatrix::execution_space >::team_size_recommended(op,vector_length);
#else
  const int team_size = Kokkos::TeamPolicy<typename AMatrix::execution_space>(rows_per_thread, Kokkos::AUTO, vector_length).team_size_recommended(op, Kokkos::ParallelForTag());
#endif
    const int rows_per_team = rows_per_thread * team_size;
    const size_type nteams = (nrow+rows_per_team-1)/rows_per_team;
    Kokkos::parallel_for("KokkosSparse::spmv<MV,NoTranspose>", Kokkos::TeamPolicy< typename AMatrix::execution_space >
       ( nteams , team_size , vector_length ) , op );

#else // KOKKOS_FAST_COMPILE this will only instantiate one Kernel for alpha/beta

    typedef SPMV_MV_LayoutLeft_Functor<AMatrix, XVector, YVector,
      2, 2, conjugate> OpType;

    typename AMatrix::const_ordinal_type nrow = A.numRows();

    OpType op (alpha, A, x, beta, y, RowsPerThread<typename AMatrix::execution_space> (NNZPerRow));

    // FIXME (mfh 07 Jun 2016) Shouldn't we use ordinal_type here
    // instead of int?  For example, if the number of threads is 1,
    // then this is just the number of rows.  Ditto for rows_per_team.
    // team_size is a hardware resource thing so it might legitimately
    // be int.
    const int rows_per_thread = RowsPerThread<typename AMatrix::execution_space >(NNZPerRow);
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    const int team_size = Kokkos::TeamPolicy< typename AMatrix::execution_space >::team_size_recommended(op,vector_length);
#else
  const int team_size = Kokkos::TeamPolicy<typename AMatrix::execution_space>(rows_per_thread, Kokkos::AUTO, vector_length).team_size_recommended(op, Kokkos::ParallelForTag());
#endif
    const int rows_per_team = rows_per_thread * team_size;
    const size_type nteams = (nrow+rows_per_team-1)/rows_per_team;
    Kokkos::parallel_for("KokkosSparse::spmv<MV,NoTranspose>",  Kokkos::TeamPolicy< typename AMatrix::execution_space >
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
static void
spmv_alpha_beta_mv_transpose (const typename YVector::non_const_value_type& alpha,
                              const AMatrix& A,
                              const XVector& x,
                              const typename YVector::non_const_value_type& beta,
                              const YVector& y)
{
  typedef typename AMatrix::ordinal_type ordinal_type;

  if (A.numRows () <= static_cast<ordinal_type> (0)) {
    return;
  }

  // We need to scale y first ("scaling" by zero just means filling
  // with zeros), since the functor works by atomic-adding into y.
  if (dobeta != 1) {
    KokkosBlas::scal (y, beta, y);
  }

  if (doalpha != 0) {
    typedef typename AMatrix::size_type size_type;

    // Assuming that no row contains duplicate entries, NNZPerRow
    // cannot be more than the number of columns of the matrix.  Thus,
    // the appropriate type is ordinal_type.
    const ordinal_type NNZPerRow = static_cast<ordinal_type> (A.nnz () / A.numRows ());

    int vector_length = 1;
    while( (static_cast<ordinal_type> (vector_length*2*3) <= NNZPerRow) && (vector_length<8) ) vector_length*=2;

#ifndef KOKKOS_FAST_COMPILE // This uses templated functions on doalpha and dobeta and will produce 16 kernels

    typedef SPMV_MV_Transpose_Functor<AMatrix, XVector, YVector,
      doalpha, dobeta, conjugate> OpType;
    OpType op (alpha, A, x, beta, y, RowsPerThread<typename AMatrix::execution_space> (NNZPerRow));

    typename AMatrix::const_ordinal_type nrow = A.numRows();

    // FIXME (mfh 07 Jun 2016) Shouldn't we use ordinal_type here
    // instead of int?  For example, if the number of threads is 1,
    // then this is just the number of rows.  Ditto for rows_per_team.
    // team_size is a hardware resource thing so it might legitimately
    // be int.
    const int rows_per_thread = RowsPerThread<typename AMatrix::execution_space >(NNZPerRow);
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    const int team_size = Kokkos::TeamPolicy< typename AMatrix::execution_space >::team_size_recommended(op,vector_length);
#else
  const int team_size = Kokkos::TeamPolicy<typename AMatrix::execution_space>(rows_per_thread, Kokkos::AUTO, vector_length).team_size_recommended(op, Kokkos::ParallelForTag());
#endif
    const int rows_per_team = rows_per_thread * team_size;
    const size_type nteams = (nrow+rows_per_team-1)/rows_per_team;
    Kokkos::parallel_for ("KokkosSparse::spmv<MV,Transpose>",  Kokkos::TeamPolicy< typename AMatrix::execution_space >
       ( nteams , team_size , vector_length ) , op );

#else // KOKKOS_FAST_COMPILE this will only instantiate one Kernel for alpha/beta

    typedef SPMV_MV_Transpose_Functor<AMatrix, XVector, YVector,
      2, 2, conjugate, SizeType> OpType;

    typename AMatrix::const_ordinal_type nrow = A.numRows();

    OpType op (alpha, A, x, beta, y, RowsPerThread<typename AMatrix::execution_space> (NNZPerRow));

    // FIXME (mfh 07 Jun 2016) Shouldn't we use ordinal_type here
    // instead of int?  For example, if the number of threads is 1,
    // then this is just the number of rows.  Ditto for rows_per_team.
    // team_size is a hardware resource thing so it might legitimately
    // be int.
    const int rows_per_thread = RowsPerThread<typename AMatrix::execution_space >(NNZPerRow);
#ifdef KOKKOS_ENABLE_DEPRECATED_CODE
    const int team_size = Kokkos::TeamPolicy< typename AMatrix::execution_space >::team_size_recommended(op,vector_length);
#else
  const int team_size = Kokkos::TeamPolicy<typename AMatrix::execution_space>(rows_per_thread, Kokkos::AUTO, vector_length).team_size_recommended(op, Kokkos::ParallelForTag());
#endif
    const int rows_per_team = rows_per_thread * team_size;
    const size_type nteams = (nrow+rows_per_team-1)/rows_per_team;
    Kokkos::parallel_for("KokkosSparse::spmv<MV,Transpose>",  Kokkos::TeamPolicy< typename AMatrix::execution_space >
       ( nteams , team_size , vector_length ) , op );

#endif // KOKKOS_FAST_COMPILE
  }
}

template<class AMatrix,
         class XVector,
         class YVector,
         int doalpha,
         int dobeta>
static void
spmv_alpha_beta_mv (const char mode[],
                    const typename YVector::non_const_value_type& alpha,
                    const AMatrix& A,
                    const XVector& x,
                    const typename YVector::non_const_value_type& beta,
                    const YVector& y)
{
  if (mode[0] == NoTranspose[0]) {
    spmv_alpha_beta_mv_no_transpose<AMatrix, XVector, YVector, doalpha, dobeta, false> (alpha, A, x, beta, y);
  }
  else if (mode[0] == Conjugate[0]) {
    spmv_alpha_beta_mv_no_transpose<AMatrix, XVector, YVector, doalpha, dobeta, true> (alpha, A, x, beta, y);
  }
  else if (mode[0] == Transpose[0]) {
    spmv_alpha_beta_mv_transpose<AMatrix, XVector, YVector, doalpha, dobeta, false> (alpha, A, x, beta, y);
  }
  else if (mode[0] == ConjugateTranspose[0]) {
    spmv_alpha_beta_mv_transpose<AMatrix, XVector, YVector, doalpha, dobeta, true> (alpha, A, x, beta, y);
  }
  else {
    Kokkos::Impl::throw_runtime_exception ("Invalid Transpose Mode for KokkosSparse::spmv()");
  }
}

template<class AMatrix,
         class XVector,
         class YVector,
         int doalpha>
void
spmv_alpha_mv (const char mode[],
               const typename YVector::non_const_value_type& alpha,
               const AMatrix& A,
               const XVector& x,
               const typename YVector::non_const_value_type& beta,
               const YVector& y)
{
  typedef typename YVector::non_const_value_type coefficient_type;
  typedef Kokkos::Details::ArithTraits<coefficient_type> KAT;

  if (beta == KAT::zero ()) {
    spmv_alpha_beta_mv<AMatrix, XVector, YVector, doalpha, 0> (mode, alpha, A, x, beta, y);
  }
  else if (beta == KAT::one ()) {
    spmv_alpha_beta_mv<AMatrix, XVector, YVector, doalpha, 1> (mode, alpha, A, x, beta, y);
  }
  else if (beta == -KAT::one ()) {
    spmv_alpha_beta_mv<AMatrix, XVector, YVector, doalpha, -1> (mode, alpha, A, x, beta, y);
  }
  else {
    spmv_alpha_beta_mv<AMatrix, XVector, YVector, doalpha, 2> (mode, alpha, A, x, beta, y);
  }
}

}
}

#endif // KOKKOSSPARSE_IMPL_SPMV_DEF_HPP_
