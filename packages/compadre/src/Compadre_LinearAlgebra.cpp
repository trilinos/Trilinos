// @HEADER
// *****************************************************************************
//     Compadre: COMpatible PArticle Discretization and REmap Toolkit
//
// Copyright 2018 NTESS and the Compadre contributors.
// SPDX-License-Identifier: BSD-2-Clause
// *****************************************************************************
// @HEADER
#include "Compadre_LinearAlgebra_Definitions.hpp"

#include "KokkosBatched_UTV_Decl.hpp"
#include "KokkosBatched_Trsm_Decl.hpp"
#include "KokkosBatched_SolveUTV_Decl_Compadre.hpp"
#include "KokkosBatched_Trsm_TeamVector_Impl.hpp"
#include "KokkosKernels_config.h"
#include <cstdlib>

// Kokkos Kernels 5.0 moved ArithTraits from namespace Kokkos (Kokkos_ArithTraits.hpp)
// to namespace KokkosKernels (KokkosKernels_ArithTraits.hpp). The old header only compiles
// with KOKKOS_ENABLE_DEPRECATED_CODE_5.
#if KOKKOSKERNELS_VERSION_MAJOR >= 5
#include "KokkosKernels_ArithTraits.hpp"
template <typename T>
using CompadreArithTraits = KokkosKernels::ArithTraits<T>;
#else
#include "Kokkos_ArithTraits.hpp"
template <typename T>
using CompadreArithTraits = Kokkos::ArithTraits<T>;
#endif

using namespace KokkosBatched;


template<class MemberType, class AViewType>
KOKKOS_INLINE_FUNCTION
int team_chol_lower_inplace(const MemberType& member, const AViewType& A, const int n) {
  using value_type = typename AViewType::non_const_value_type;
  // a pivot at or below this fraction of its original diagonal entry means A is not positive definite
  const value_type relative_pivot_tolerance = 10 * n * CompadreArithTraits<value_type>::epsilon();
  for (int j=0; j<n; ++j) {

    // A(j,j) -= sum_{k<j} A(j,k)^2
    typename AViewType::non_const_value_type s = 0;
    Kokkos::parallel_reduce(Kokkos::ThreadVectorRange(member, j), [&](int k, auto& lsum) {
      lsum += A(j,k)*A(j,k);
    }, s);
    member.team_barrier();

    if (member.team_rank() == 0) {
      const value_type a_jj = A(j,j);
      const value_type pivot = a_jj - s;
      // A = (P*sqrt(W))^T*(P*sqrt(W)) is positive definite only if P*sqrt(W) has full column rank.
      // A pivot that vanishes relative to its original diagonal entry (or is negative or NaN, for which
      // the comparison is false) means the problem is rank deficient.
      if (pivot > relative_pivot_tolerance * a_jj) {
        A(j,j) = Kokkos::sqrt(pivot);
      } else {
        Kokkos::abort("Compadre: CHOL solver failed because P^T*W*P is not positive definite, i.e. the "
                      "GMLS problem is rank deficient for at least one target. Either:\n"
                      "  a.) provide more neighbors for each target, or\n"
                      "  b.) change the sampling operator, as it may have a null space over the reconstruction\n"
                      "      space (e.g. sampling a gradient cannot recover a constant), or\n"
                      "  c.) use the QR solver, which handles rank deficient problems.\n");
      }
    }
    member.team_barrier();

    // For i=j+1..n-1: A(i,j) = (A(i,j) - sum_{k<j} A(i,k)A(j,k))/A(j,j)
    Kokkos::parallel_for(Kokkos::TeamThreadRange(member, j+1, n), [&](int i) {
      typename AViewType::non_const_value_type t = 0;
      for (int k=0; k<j; ++k) t += A(i,k)*A(j,k);
      A(i,j) = (A(i,j) - t)/A(j,j);
    });
    member.team_barrier();
  }

  // we do not zero the upper triangle
  return 0;
}

namespace Compadre{
namespace GMLS_LinearAlgebra {

  template<typename DeviceType,
           typename AlgoTagType,
           typename MatrixViewType_A,
           typename MatrixViewType_B,
           typename MatrixViewType_X>
  struct Functor_TestBatchedTeamVectorSolveUTV {
    MatrixViewType_A _a;
    MatrixViewType_B _b;

    int _pm_getTeamScratchLevel_0;
    int _pm_getTeamScratchLevel_1;
    int _M, _N, _NRHS;
    bool _implicit_RHS;

    KOKKOS_INLINE_FUNCTION
    Functor_TestBatchedTeamVectorSolveUTV(
                      const int M,
                      const int N,
                      const int NRHS,
                      const MatrixViewType_A &a,
                      const MatrixViewType_B &b,
                      const bool implicit_RHS)
      : _a(a), _b(b), _M(M), _N(N), _NRHS(NRHS), _implicit_RHS(implicit_RHS) 
        { _pm_getTeamScratchLevel_0 = 0; _pm_getTeamScratchLevel_1 = 0; }

    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void operator()(const MemberType &member) const {

      const int k = member.league_rank();

      // workspace vectors
      scratch_vector_type ww_fast(member.team_scratch(_pm_getTeamScratchLevel_0), 3*_M);
      scratch_vector_type ww_slow(member.team_scratch(_pm_getTeamScratchLevel_1), TO_GLOBAL(_N)*TO_GLOBAL(_NRHS));

      device_unmanaged_matrix_right_type aa(_a.data() + TO_GLOBAL(k)*TO_GLOBAL(_a.extent(1))*TO_GLOBAL(_a.extent(2)), 
              _a.extent(1), _a.extent(2));
      device_unmanaged_matrix_right_type bb(_b.data() + TO_GLOBAL(k)*TO_GLOBAL(_b.extent(1))*TO_GLOBAL(_b.extent(2)), 
              _b.extent(1), _b.extent(2));
      device_unmanaged_matrix_right_type xx(_b.data() + TO_GLOBAL(k)*TO_GLOBAL(_b.extent(1))*TO_GLOBAL(_b.extent(2)), 
              _b.extent(1), _b.extent(2));

      // if sizes don't match extents, then copy to a view with extents matching sizes
      if ((size_t)_M!=_a.extent(1) || (size_t)_N!=_a.extent(2)) {
        scratch_matrix_right_type tmp(ww_slow.data(), _M, _N);
        auto aaa = device_unmanaged_matrix_right_type(_a.data() + TO_GLOBAL(k)*TO_GLOBAL(_a.extent(1))*TO_GLOBAL(_a.extent(2)), _M, _N);
        // copy A to W, then back to A
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,_M),[&](const int &i) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,0,_N),[&](const int &j) {
              tmp(i,j) = aa(i,j);
          });
        });
        member.team_barrier();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,_M),[&](const int &i) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,0,_N),[&](const int &j) {
              aaa(i,j) = tmp(i,j);
          });
        });
        member.team_barrier();
        aa = aaa;
      }

      if (std::is_same<typename MatrixViewType_B::array_layout, layout_left>::value) {
        scratch_matrix_right_type tmp(ww_slow.data(), _N, _NRHS);
        // coming from CHOL
        // then copy B to W, then back to B
        auto bb_left = 
            device_unmanaged_matrix_left_type(_b.data() + TO_GLOBAL(k)*TO_GLOBAL(_b.extent(1))*TO_GLOBAL(_b.extent(2)), 
                    _b.extent(1), _b.extent(2));
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,_N),[&](const int &i) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,0,_NRHS),[&](const int &j) {
              tmp(i,j) = bb_left(i,j);
          });
        });
        member.team_barrier();
        Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,_N),[&](const int &i) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,0,_NRHS),[&](const int &j) {
              bb(i,j) = tmp(i,j);
          });
        });
      }

      scratch_matrix_right_type uu(member.team_scratch(_pm_getTeamScratchLevel_1), _M, _N /* only N columns of U are filled, maximum */);
      scratch_matrix_right_type vv(member.team_scratch(_pm_getTeamScratchLevel_1), _N, _N);
      scratch_local_index_type pp(member.team_scratch(_pm_getTeamScratchLevel_0), _N);

      bool do_print = false;
      if (do_print) {
        Kokkos::single(Kokkos::PerTeam(member), [&] () {
          using Kokkos::printf;
          //print a
          printf("a=zeros(%lu,%lu);\n", aa.extent(0), aa.extent(1));
              for (size_t i=0; i<aa.extent(0); ++i) {
                  for (size_t j=0; j<aa.extent(1); ++j) {
                      printf("a(%lu,%lu)= %f;\n", i+1,j+1, aa(i,j));
                  }
              }
          //print b
          printf("b=zeros(%lu,%lu);\n", bb.extent(0), bb.extent(1));
              for (size_t i=0; i<bb.extent(0); ++i) {
                  for (size_t j=0; j<bb.extent(1); ++j) {
                      printf("b(%lu,%lu)= %f;\n", i+1,j+1, bb(i,j));
                  }
              }
        });
      }
      do_print = false;

      /// Solving Ax = b using UTV transformation
      /// A P^T P x = b
      /// UTV P x = b;

      /// UTV = A P^T
      int matrix_rank(0);
      member.team_barrier();
      TeamVectorUTV<MemberType,AlgoTagType>
        ::invoke(member, aa, pp, uu, vv, ww_fast, matrix_rank);
      member.team_barrier();

      if (do_print) {
        Kokkos::single(Kokkos::PerTeam(member), [&] () {
        using Kokkos::printf;
        printf("matrix_rank: %d\n", matrix_rank);
        //print u
        printf("u=zeros(%lu,%lu);\n", uu.extent(0), uu.extent(1));
        for (size_t i=0; i<uu.extent(0); ++i) {
            for (size_t j=0; j<uu.extent(1); ++j) {
                printf("u(%lu,%lu)= %f;\n", i+1,j+1, uu(i,j));
            }
        }
        });
      }

      TeamVectorSolveUTVCompadre<MemberType,AlgoTagType>
        ::invoke(member, matrix_rank, _M, _N, _NRHS, uu, aa, vv, pp, bb, xx, ww_slow, ww_fast, _implicit_RHS);
      member.team_barrier();

    }

    inline
    void run(ParallelManager pm) {
      typedef typename MatrixViewType_A::non_const_value_type value_type;
      std::string name_region("KokkosBatched::Test::TeamVectorSolveUTVCompadre");
      std::string name_value_type = ( std::is_same<value_type,float>::value ? "::Float" :
                                      std::is_same<value_type,double>::value ? "::Double" :
                                      std::is_same<value_type,Kokkos::complex<float> >::value ? "::ComplexFloat" :
                                      std::is_same<value_type,Kokkos::complex<double> >::value ? "::ComplexDouble" : "::UnknownValueType" );
      std::string name = name_region + name_value_type;
      Kokkos::Profiling::pushRegion( name.c_str() );

      _pm_getTeamScratchLevel_0 = pm.getTeamScratchLevel(0);
      _pm_getTeamScratchLevel_1 = pm.getTeamScratchLevel(1);
      
      int scratch_size = scratch_matrix_right_type::shmem_size(_N, _N); // V
      scratch_size += scratch_matrix_right_type::shmem_size(_M, _N /* only N columns of U are filled, maximum */); // U
      scratch_size += scratch_vector_type::shmem_size(TO_GLOBAL(_N)*TO_GLOBAL(_NRHS)); // W (for SolveUTV)

      int l0_scratch_size = scratch_vector_type::shmem_size(_N); // P (temporary)
      l0_scratch_size += scratch_vector_type::shmem_size(3*_M); // W (for UTV)

      pm.clearScratchSizes();
      pm.setTeamScratchSize(0, l0_scratch_size);
      pm.setTeamScratchSize(1, scratch_size);

      auto tp = pm.TeamPolicyThreadsAndVectors(_a.extent(0));
      // optionally solve each matrix serially (1 thread, 1 vector lane) when N < $QR_SERIAL_BELOW_N
      if (const char* qr_serial_below_n = std::getenv("QR_SERIAL_BELOW_N")) {
        if (_N < std::atoi(qr_serial_below_n)) {
          tp = pm.TeamPolicyThreadsAndVectors(_a.extent(0), 1, 1);
        }
      }
      Kokkos::parallel_for("BatchedTeamVectorSolveUTV", tp, *this);
      Kokkos::fence();

      Kokkos::Profiling::popRegion();
    }
  };



  template<typename DeviceType,
           typename AlgoTagType,
           typename MatrixViewType_A,
           typename MatrixViewType_B,
           typename MatrixViewType_X>
  struct Functor_TestBatchedTeamVectorSolveCholesky {
    MatrixViewType_A _a;
    MatrixViewType_B _b;

    int _M, _N, _NRHS;

    KOKKOS_INLINE_FUNCTION
    Functor_TestBatchedTeamVectorSolveCholesky(
                      const int M,
                      const int N,
                      const int NRHS,
                      const MatrixViewType_A &a,
                      const MatrixViewType_B &b)
      : _a(a), _b(b), _M(M), _N(N), _NRHS(NRHS) {}

    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void operator()(const MemberType &member) const {

      const int k = member.league_rank();

      device_unmanaged_matrix_right_type aa(_a.data() + TO_GLOBAL(k)*TO_GLOBAL(_a.extent(1))*TO_GLOBAL(_a.extent(2)), 
              _a.extent(1), _a.extent(2));
      device_unmanaged_matrix_right_type bb(_b.data() + TO_GLOBAL(k)*TO_GLOBAL(_b.extent(1))*TO_GLOBAL(_b.extent(2)), 
              _b.extent(1), _b.extent(2));
      device_unmanaged_matrix_right_type xx(_b.data() + TO_GLOBAL(k)*TO_GLOBAL(_b.extent(1))*TO_GLOBAL(_b.extent(2)), 
              _b.extent(1), _b.extent(2));

      auto a_small = Kokkos::subview(aa,
                              Kokkos::make_pair<size_t,size_t>(0, _M),
                              Kokkos::make_pair<size_t,size_t>(0, _N));

      auto x_small = Kokkos::subview(xx,
                              Kokkos::make_pair<size_t,size_t>(0, _N),
                              Kokkos::make_pair<size_t,size_t>(0, _NRHS));

      team_chol_lower_inplace(member, a_small, _N);
      member.team_barrier();

      using Trsmty = KokkosBatched::TeamVectorTrsm<
          MemberType,
          KokkosBatched::Side::Left,
          KokkosBatched::Uplo::Lower,
          KokkosBatched::Trans::NoTranspose,
          KokkosBatched::Diag::NonUnit,
          AlgoTagType>;
      
      Trsmty::invoke(member, 1.0, a_small, x_small);
      member.team_barrier();
      
      using TrsmT = KokkosBatched::TeamVectorTrsm<
          MemberType,
          KokkosBatched::Side::Left,
          KokkosBatched::Uplo::Lower,
          KokkosBatched::Trans::Transpose,
          KokkosBatched::Diag::NonUnit,
          AlgoTagType>;
      
      TrsmT::invoke(member, 1.0, a_small, x_small);
      member.team_barrier();

    }

    inline
    void run(ParallelManager pm) {
      typedef typename MatrixViewType_A::non_const_value_type value_type;
      std::string name_region("KokkosBatched::Test::TeamVectorSolveCholesky");
      std::string name_value_type = ( std::is_same<value_type,float>::value ? "::Float" :
                                      std::is_same<value_type,double>::value ? "::Double" :
                                      std::is_same<value_type,Kokkos::complex<float> >::value ? "::ComplexFloat" :
                                      std::is_same<value_type,Kokkos::complex<double> >::value ? "::ComplexDouble" : "::UnknownValueType" );
      std::string name = name_region + name_value_type;
      Kokkos::Profiling::pushRegion( name.c_str() );

      pm.clearScratchSizes();
      auto tp = pm.TeamPolicyThreadsAndVectors(_a.extent(0));
      // optionally solve each matrix serially (1 thread, 1 vector lane) when N < $CHOL_SERIAL_BELOW_N
      if (const char* chol_serial_below_n = std::getenv("CHOL_SERIAL_BELOW_N")) {
        if (_N < std::atoi(chol_serial_below_n)) {
          tp = pm.TeamPolicyThreadsAndVectors(_a.extent(0), 1, 1);
        }
      }
      Kokkos::parallel_for("BatchedTeamVectorSolveCholesky", tp, *this);
      Kokkos::fence();

      Kokkos::Profiling::popRegion();
    }
  };



  /*! \brief Transpose a matrix of original shape N x M (padded). Result will be M x N.
      \param M       [in] - number of rows in the result
      \param N       [in] - number of columns in the result
      \param b   [in/out] - in as rank 3 matrix (num_matrices * N * M), out as rank 3 matrix (num_matrices * M * N)
  */
  template<typename DeviceType,
           typename MatrixViewType>
  struct Functor_Transpose {
    MatrixViewType _b;
    int _M, _N;
    int _pm_getTeamScratchLevel_1;

    KOKKOS_INLINE_FUNCTION
    Functor_Transpose(
                      const int M,
                      const int N,
                      const MatrixViewType &b)
      : _b(b), _M(M), _N(N) {}

    template<typename MemberType>
    KOKKOS_INLINE_FUNCTION
    void operator()(const MemberType &member) const {

      const int k = member.league_rank();

      // workspace vector
      scratch_vector_type ww_slow(member.team_scratch(_pm_getTeamScratchLevel_1), TO_GLOBAL(_M)*TO_GLOBAL(_N));
      device_unmanaged_matrix_right_type bb(_b.data() + TO_GLOBAL(k)*TO_GLOBAL(_b.extent(1))*TO_GLOBAL(_b.extent(2)), 
              _b.extent(1), _b.extent(2));

      scratch_matrix_right_type tmp(ww_slow.data(), _M, _N);
      // coming from CHOL
      // then copy B to W, then back to B
      auto bb_left = 
          device_unmanaged_matrix_left_type(_b.data() + TO_GLOBAL(k)*TO_GLOBAL(_b.extent(1))*TO_GLOBAL(_b.extent(2)), 
                  _b.extent(1), _b.extent(2));
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,_M),[&](const int &i) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,0,_N),[&](const int &j) {
            tmp(i,j) = bb_left(i,j);
        });
      });
      member.team_barrier();
      Kokkos::parallel_for(Kokkos::TeamThreadRange(member,0,_M),[&](const int &i) {
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member,0,_N),[&](const int &j) {
            bb(i,j) = tmp(i,j);
        });
      });
      member.team_barrier();
    }

    inline
    void run(ParallelManager pm) {
      typedef typename MatrixViewType::non_const_value_type value_type;
      std::string name_region("KokkosBatched::Test::Transpose");
      std::string name_value_type = ( std::is_same<value_type,float>::value ? "::Float" :
                                      std::is_same<value_type,double>::value ? "::Double" :
                                      std::is_same<value_type,Kokkos::complex<float> >::value ? "::ComplexFloat" :
                                      std::is_same<value_type,Kokkos::complex<double> >::value ? "::ComplexDouble" : "::UnknownValueType" );
      std::string name = name_region + name_value_type;
      Kokkos::Profiling::pushRegion( name.c_str() );

      _pm_getTeamScratchLevel_1 = pm.getTeamScratchLevel(1);
      
      size_t scratch_size = scratch_vector_type::shmem_size(TO_GLOBAL(_M)*TO_GLOBAL(_N));

      pm.clearScratchSizes();
      pm.setTeamScratchSize(1, scratch_size);

      // set up number of threads and vectors based on size
      auto tp = pm.TeamPolicyThreadsAndVectors(_b.extent(0));
      Kokkos::parallel_for("BatchedTeamVectorTranspose", tp, *this);
      Kokkos::fence();

      Kokkos::Profiling::popRegion();
    }
  };

void batchSolve(ParallelManager pm, DenseSolverType dst, double *A, int lda, int nda, bool transpose_A, double *B, int ldb, int ndb, bool transpose_B, int M, int N, int NRHS, const int num_matrices, const bool implicit_RHS, const bool rank_full) {

    typedef Algo::UTV::Unblocked algo_tag_type;
    typedef Kokkos::View<double***, layout_right, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
                    MatrixViewType_A;
    typedef Kokkos::View<double***, layout_right, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
                    MatrixViewType_B;
    typedef Kokkos::View<double***, layout_right, Kokkos::MemoryTraits<Kokkos::Unmanaged> >
                    MatrixViewType_X;


    MatrixViewType_A mat_A(A, num_matrices, lda, nda);
    MatrixViewType_B mat_B(B, num_matrices, ldb, ndb);

    // an implicit RHS stores the diagonal of sqrt(W) as a vector, so it cannot be transposed
    compadre_assert_debug((!(transpose_B && implicit_RHS)) && "transpose_B is not supported when implicit_RHS is true.");

    // important: A and B are rank 3, with potential padding around their 2nd and 3rd rank
    // so it is not necessarily tru that lda or nda == M, for instance
    //
    // sizes given to tranpose are w.r.t. the desired, post tranpose layout
    if (transpose_A) {
        // so e.g., N is related to lda and M is related to nda
        mat_A = MatrixViewType_A(A, num_matrices, nda, lda);
        Functor_Transpose
          <device_execution_space, MatrixViewType_A>(M,N,mat_A).run(pm);
    }
    if (transpose_B) {
        // so e.g., N is related to lda and M is related to nda
        // B is P, but we need P^T
        mat_B = MatrixViewType_B(B, num_matrices, ndb, ldb);
        Functor_Transpose
          <device_execution_space, MatrixViewType_B>(N,NRHS,mat_B).run(pm);
    }

    if (dst==DenseSolverType::CHOL && rank_full) {
        Functor_TestBatchedTeamVectorSolveCholesky
          <device_execution_space, algo_tag_type, MatrixViewType_A, MatrixViewType_B, MatrixViewType_X>(M,N,NRHS,mat_A,mat_B).run(pm);
    } else {
        Functor_TestBatchedTeamVectorSolveUTV
          <device_execution_space, algo_tag_type, MatrixViewType_A, MatrixViewType_B, MatrixViewType_X>(M,N,NRHS,mat_A,mat_B,implicit_RHS).run(pm);
    }
    Kokkos::fence();

}


} // GMLS_LinearAlgebra
} // Compadre
