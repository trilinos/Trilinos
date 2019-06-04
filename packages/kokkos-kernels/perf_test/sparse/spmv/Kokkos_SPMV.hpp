/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef KOKKOS_SPMV_HPP_
#define KOKKOS_SPMV_HPP_



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

  KOKKOS_INLINE_FUNCTION void
  operator() (const team_member& dev) const
  {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(dev,0,rows_per_team), [&] (const ordinal_type& loop) {

      const ordinal_type iRow = static_cast<ordinal_type> ( dev.league_rank() ) * rows_per_team + loop;
      if (iRow >= m_A.numRows ()) {
        return;
      }
      const KokkosSparse::SparseRowViewConst<AMatrix> row = m_A.rowConst(iRow);
      const ordinal_type row_length = static_cast<ordinal_type> (row.length);
      value_type sum = 0;

      Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(dev,row_length), [&] (const ordinal_type& iEntry, value_type& lsum) {
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
int launch_parameters(int numRows, int nnz, int rows_per_thread, int& team_size, int& vector_length) {
  int rows_per_team;
  int nnz_per_row = nnz/numRows;
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
  if(team_size < 1)
    team_size = 256/vector_length;
  #endif

  rows_per_team = rows_per_thread * team_size;

  if(rows_per_team < 0) {
    int nnz_per_team = 4096;
    int conc = execution_space::concurrency();
    while((conc * nnz_per_team * 4> nnz)&&(nnz_per_team>256)) nnz_per_team/=2;
    int tmp_nnz_per_row = nnz/numRows;
    rows_per_team = (nnz_per_team+tmp_nnz_per_row - 1)/tmp_nnz_per_row;
  }


  return rows_per_team;
}

template<typename AType, typename XType, typename YType,class ScheduleType>
void kk_matvec(AType A, XType x, YType y, int rows_per_thread, int team_size, int vector_length) {

  typedef typename XType::non_const_value_type Scalar;
  typedef typename AType::execution_space execution_space;
  typedef KokkosSparse::CrsMatrix<const Scalar,int,execution_space,void,int> matrix_type ;
  typedef typename Kokkos::View<Scalar*,Kokkos::LayoutLeft,execution_space> y_type;
  typedef typename Kokkos::View<const Scalar*,Kokkos::LayoutLeft,execution_space,Kokkos::MemoryRandomAccess > x_type;

  int rows_per_team = launch_parameters<execution_space>(A.numRows(),A.nnz(),rows_per_thread,team_size,vector_length);

  double s_a = 1.0;
  double s_b = 0.0;
  SPMV_Functor<matrix_type,x_type,y_type,0,false> func (s_a,A,x,s_b,y,rows_per_team);

  int worksets = (y.extent(0)+rows_per_team-1)/rows_per_team;

  Kokkos::TeamPolicy<Kokkos::Schedule<ScheduleType> > policy(1,1);

  if(team_size>0)
    policy = Kokkos::TeamPolicy<Kokkos::Schedule<ScheduleType> >(worksets,team_size,vector_length);
  else
    policy = Kokkos::TeamPolicy<Kokkos::Schedule<ScheduleType> >(worksets,Kokkos::AUTO,vector_length);

  Kokkos::parallel_for("KokkosSparse::PerfTest::SpMV", policy,func);
}


#endif /* KOKKOS_SPMV_HPP_ */
