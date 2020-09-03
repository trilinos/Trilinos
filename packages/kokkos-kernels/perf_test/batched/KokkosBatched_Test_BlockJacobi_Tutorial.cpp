/// Kokkos headers
#include "Kokkos_Core.hpp"
#include "impl/Kokkos_Timer.hpp"
#include "Kokkos_Random.hpp"


#if  defined(KOKKOS_ENABLE_CXX11_DISPATCH_LAMBDA)
#if !defined(KOKKOS_ENABLE_CUDA) || (8000 <= CUDA_VERSION)
#if  defined(KOKKOS_ENABLE_CUDA_LAMBDA)
#define KOKKOSBATCHED_TEST_BLOCKJACOBI 
#endif 
#endif
#endif

#if defined(KOKKOSBATCHED_TEST_BLOCKJACOBI) 

/// KokkosKernels headers
#include "KokkosBatched_Util.hpp"

#include "Kokkos_ArithTraits.hpp"
#include "KokkosBatched_Util.hpp"
#include "KokkosBatched_Copy_Decl.hpp"
#include "KokkosBatched_Copy_Impl.hpp"
#include "KokkosBatched_SetIdentity_Decl.hpp"
#include "KokkosBatched_SetIdentity_Impl.hpp"
#include "KokkosBatched_Gemv_Decl.hpp"
#include "KokkosBatched_Gemv_Team_Impl.hpp"
#include "KokkosBatched_Trsm_Decl.hpp"
#include "KokkosBatched_Trsm_Team_Impl.hpp"
#include "KokkosBatched_LU_Decl.hpp"
#include "KokkosBatched_LU_Team_Impl.hpp"

/// cuda profile
#if defined(KOKKOS_ENABLE_CUDA)
#include "cuda_profiler_api.h"
#endif


using exec_space_type =  Kokkos::DefaultExecutionSpace;
using memory_space_type = typename exec_space_type::memory_space;
using host_space = Kokkos::DefaultHostExecutionSpace;

using val_type = double;
using policy_type = Kokkos::TeamPolicy<exec_space_type>;
using member_type = typename policy_type::member_type;

using namespace KokkosBatched;

template<typename ManyMatrixType,
	 typename ManyVectorType>
val_type computeResidual(const ManyMatrixType &A,
			   const ManyVectorType &x,
			   const ManyVectorType &b,
			   const ManyVectorType &r) {
  /// compute residual
  val_type residual(0);
  {
    policy_type policy(A.extent(0), Kokkos::AUTO());	
    Kokkos::deep_copy(r, b);
    Kokkos::parallel_reduce
      ("compute-residual",
       policy, KOKKOS_LAMBDA(const member_type &member, val_type &update) {
	const int i = member.league_rank();
	const val_type one(1);
	auto AA = Kokkos::subview(A, i, Kokkos::ALL(), Kokkos::ALL());
	auto xx = Kokkos::subview(x, i, Kokkos::ALL());
	auto rr = Kokkos::subview(r, i, Kokkos::ALL());
	
	TeamGemv<member_type,
		 Trans::NoTranspose,
		 Algo::Level2::Unblocked>
	  ::invoke(member, -one, AA, xx, one, rr);

        val_type sum(0);
	Kokkos::parallel_reduce
	  (Kokkos::TeamThreadRange(member, rr.extent(0)),
	   [&](const int &k, val_type &lsum) {
            lsum += Kokkos::ArithTraits<val_type>::abs(rr(k));
          }, sum);
        Kokkos::single(Kokkos::PerTeam(member), [&]() {
            update += sum;
          });
      }, residual);
  }
  return residual;
}

int main(int argc, char* argv[]) {
  Kokkos::initialize(argc, argv);
  {
#if defined(KOKKOS_ENABLE_CUDA)
    cudaProfilerStop();
#endif
    Kokkos::print_configuration(std::cout);

    Kokkos::Impl::Timer timer;

    ///
    /// input arguments parsing
    ///
    int N = 128*128; /// # of problems (batch size)
    int Blk = 5;     /// block dimension
    for (int i=1;i<argc;++i) {
      const std::string& token = argv[i];
      if (token == std::string("-N")) N = std::atoi(argv[++i]);
      if (token == std::string("-B")) Blk = std::atoi(argv[++i]);
    }
    printf(" :::: Testing (N = %d, Blk = %d)\n", N, Blk);


    ///
    /// Problem container: rank-3 array
    ///
    /// A - multiple block matrices representing block diagonals
    /// T - temporal block matrices to store its LU factors
    /// x - solution vector
    /// b - right hand side vector
    ///
    Kokkos::View<val_type***,Kokkos::LayoutRight,exec_space_type> A("block diagonals", N, Blk, Blk);
    Kokkos::View<val_type***,Kokkos::LayoutRight,exec_space_type> T("temporal block diagonals", N, Blk, Blk);    
    Kokkos::View<val_type**,Kokkos::LayoutRight,exec_space_type> x("x", N, Blk);
    Kokkos::View<val_type**,Kokkos::LayoutRight,exec_space_type> b("b", N, Blk);

    /// copy of A to check residual
    Kokkos::View<val_type***,Kokkos::LayoutRight,exec_space_type> Acopy("Acopy",
								     A.extent(0),
								     A.extent(1),
								     A.extent(2));

    /// residual vector
    Kokkos::View<val_type**,Kokkos::LayoutRight,exec_space_type> r("r",
								b.extent(0), 
								b.extent(1));

    /// The block diagonal matrices are assumed to be extracted from a block sparse matrix.
    /// Here we set the blocks with random values
    Kokkos::Random_XorShift64_Pool<exec_space_type> random(13245);
    Kokkos::fill_random(A, random, val_type(1.0));
    Kokkos::fill_random(b, random, val_type(1.0));

    Kokkos::deep_copy(Acopy, A);

    ///
    /// Objective :
    ///    - Construct the inverse of A(i,:,:) for all i.
    ///    - Solve the equation using matrix vector multiplication. 
    
    /// Task 1. Use the so-called standard batch interface
    ///    parallel_for(factorize)
    ///    parallel_For(set identity matrix)
    ///    parallel_for(solve lower triangular)
    ///    parallel_for(solve upper triangular)
    ///    parallel_for(matrix vector multiplication)

    {
#if defined(KOKKOS_ENABLE_CUDA)
      cudaProfilerStart();
#endif
      Kokkos::deep_copy(A, Acopy);
      
      /// construction of block jacobi using batched blas interface
      /// each parallel for is a batch function
      {
	policy_type policy(A.extent(0), Kokkos::AUTO());	
	timer.reset();
	Kokkos::parallel_for
	  ("task1.factorize",
	   policy, KOKKOS_LAMBDA(const member_type &member) {
	    const int i = member.league_rank();
	    auto AA = Kokkos::subview(A, i, Kokkos::ALL(), Kokkos::ALL());
	    TeamLU<member_type,Algo::Level3::Unblocked>::invoke(member,AA);
	  });
	Kokkos::deep_copy(T, A);
	Kokkos::parallel_for
	  ("task1.set-identity",
	   policy, KOKKOS_LAMBDA(const member_type &member) {
	    const int i = member.league_rank();
	    auto AA = Kokkos::subview(A, i, Kokkos::ALL(), Kokkos::ALL());	  
	    TeamSetIdentity<member_type>::invoke(member, AA);
	  });
	Kokkos::fence();      
	Kokkos::parallel_for
	  ("task1.solve-lower-triangular",
	   policy, KOKKOS_LAMBDA(const member_type &member) {
	    const int i = member.league_rank();
	    const val_type one(1);
	    auto AA = Kokkos::subview(A, i, Kokkos::ALL(), Kokkos::ALL());
	    auto TT = Kokkos::subview(T, i, Kokkos::ALL(), Kokkos::ALL());
	    TeamTrsm<member_type,
		     Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::Unit,
		     Algo::Level3::Unblocked>
	      ::invoke(member, one, TT, AA);
	  });
	Kokkos::fence();      
	Kokkos::parallel_for
	  ("task1.solve-upper-triangular",
	   policy, KOKKOS_LAMBDA(const member_type &member) {
	    const int i = member.league_rank();
	    const val_type one(1);
	    auto AA = Kokkos::subview(A, i, Kokkos::ALL(), Kokkos::ALL());
	    auto TT = Kokkos::subview(T, i, Kokkos::ALL(), Kokkos::ALL());
	    TeamTrsm<member_type,
		     Side::Left,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,
		     Algo::Level3::Unblocked>
	      ::invoke(member, one, TT, AA);
	  });      
	Kokkos::fence();
	const double t = timer.seconds();
	printf("task 1: construction of jacobi time = %f , # of constructions per min = %.0f \n", t, 1.0/t*60);
      }
      
      /// apply block jacobi
      {
	timer.reset();
	policy_type policy(A.extent(0), Kokkos::AUTO());
	Kokkos::parallel_for
	  ("task1.apply-block-jacobi",
	   policy, KOKKOS_LAMBDA(const member_type &member) {
	    const int i = member.league_rank();
	    const val_type one(1), zero(0);
	    auto AA = Kokkos::subview(A, i, Kokkos::ALL(), Kokkos::ALL());
	    auto xx = Kokkos::subview(x, i, Kokkos::ALL());
	    auto bb = Kokkos::subview(b, i, Kokkos::ALL());
	    TeamGemv<member_type,
		     Trans::NoTranspose,
		     Algo::Level2::Unblocked>
	      ::invoke(member, one, AA, bb, zero, xx);
	  });
	const double t = timer.seconds();
	printf("task 1: application of jacobi time = %f , # of applications per min = %.0f \n", t, 1.0/t*60);	
      }

      /// check residual
      {
	const double residual = computeResidual(Acopy, x, b, r);
	printf("task 1: residual = %f\n", residual);
      }
#if defined(KOKKOS_ENABLE_CUDA) 
      cudaProfilerStop();
#endif      
    }
    
    /// Task 2. Compose a new batch function using kokkos batched team-level interface
    ///    parallel_for(LU, set identity, solve lower/upper triangular)
    ///    parallel_for(matrix vector multiplication)

    {
#if defined(KOKKOS_ENABLE_CUDA)
      cudaProfilerStart();
#endif
      Kokkos::deep_copy(A, Acopy);
      
      /// construction of block jacobi using batched blas interface
      /// each parallel for is a batch function
      {
	policy_type policy(A.extent(0), Kokkos::AUTO());	
	timer.reset();
	Kokkos::parallel_for
	  ("task2.factorize-invert",
	   policy, KOKKOS_LAMBDA(const member_type &member) {
	    const val_type one(1);
	    const int i = member.league_rank();
	    auto AA = Kokkos::subview(A, i, Kokkos::ALL(), Kokkos::ALL());
	    auto TT = Kokkos::subview(T, i, Kokkos::ALL(), Kokkos::ALL());
	    
	    TeamLU<member_type,Algo::Level3::Unblocked>::invoke(member,AA);
	    TeamCopy<member_type,Trans::NoTranspose>::invoke(member, AA, TT);
	    TeamSetIdentity<member_type>::invoke(member, AA);
	    TeamTrsm<member_type,
		     Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::Unit,
		     Algo::Level3::Unblocked>
	      ::invoke(member, one, TT, AA);
	    TeamTrsm<member_type,
		     Side::Left,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,
		     Algo::Level3::Unblocked>
	      ::invoke(member, one, TT, AA);
	  });
	Kokkos::fence();
	const double t = timer.seconds();
	printf("task 2: construction of jacobi time = %f , # of constructions per min = %.0f \n", t, 1.0/t*60);
      }
      
      /// apply block jacobi
      {
	timer.reset();
	policy_type policy(A.extent(0), Kokkos::AUTO());
	Kokkos::parallel_for
	  ("task2.apply-block-jacobi",
	   policy, KOKKOS_LAMBDA(const member_type &member) {
	    const int i = member.league_rank();
	    const val_type one(1), zero(0);
	    auto AA = Kokkos::subview(A, i, Kokkos::ALL(), Kokkos::ALL());
	    auto xx = Kokkos::subview(x, i, Kokkos::ALL());
	    auto bb = Kokkos::subview(b, i, Kokkos::ALL());
	    TeamGemv<member_type,
		     Trans::NoTranspose,
		     Algo::Level2::Unblocked>
	      ::invoke(member, one, AA, bb, zero, xx);
	  });
	const double t = timer.seconds();
	printf("task 2: application of jacobi time = %f , # of applications per min = %.0f \n", t, 1.0/t*60);	
      }

      /// check residual
      {
	const double residual = computeResidual(Acopy, x, b, r);
	printf("task 2: residual = %f\n", residual);
      }
#if defined(KOKKOS_ENABLE_CUDA) 
      cudaProfilerStop();
#endif      
    }

  }
  Kokkos::finalize();

  return 0;
}


#else
int main() {
  return 0;
}
#endif

