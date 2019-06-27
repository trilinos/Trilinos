#include <cassert>
#include <vector>
#include <algorithm>

#include "KokkosBatched_Util.hpp"
#if defined(__KOKKOSBATCHED_INTEL_MKL__)
#include "mkl.h"
#endif

#if defined(__KOKKOSBATCHED_NVIDIA_CUBLAS__)
#include "cuda_runtime.h"
#include "cublas_v2.h"
#include "cublas_api.h"
#endif

#include "Kokkos_Core.hpp"
#include "impl/Kokkos_Timer.hpp"

#include "KokkosBatched_Vector.hpp"

#include "KokkosBatched_Copy_Decl.hpp"
#include "KokkosBatched_Copy_Impl.hpp"

#include "KokkosBatched_AddRadial_Decl.hpp"
#include "KokkosBatched_AddRadial_Impl.hpp"

#include "KokkosBatched_Gemv_Decl.hpp"
#include "KokkosBatched_Gemv_Serial_Impl.hpp"
#include "KokkosBatched_Gemv_Team_Impl.hpp"

#include "KokkosBatched_Trsv_Decl.hpp"
#include "KokkosBatched_Trsv_Serial_Impl.hpp"
#include "KokkosBatched_Trsv_Team_Impl.hpp"

#include "KokkosBatched_Gemm_Decl.hpp"
#include "KokkosBatched_Gemm_Serial_Impl.hpp"
#include "KokkosBatched_Gemm_Team_Impl.hpp"

#include "KokkosBatched_Trsm_Decl.hpp"
#include "KokkosBatched_Trsm_Serial_Impl.hpp"
#include "KokkosBatched_Trsm_Team_Impl.hpp"

#include "KokkosBatched_LU_Decl.hpp"
#include "KokkosBatched_LU_Serial_Impl.hpp"
#include "KokkosBatched_LU_Team_Impl.hpp"

#include "KokkosBatched_Test_BlockCrs_Util.hpp"

namespace KokkosBatched {
  namespace Test {

    struct RangeTag {};
    struct TeamTag {};
    struct TeamShmemTag {};
    
    template <typename ExecSpace, typename ValueType, typename ArrayLayout,
              int VectorLength,
              typename LU_AlgoTagType,
              typename Trsm_AlgoTagType,
              typename Gemm_AlgoTagType>
    class FactorizeBlockTridiagMatrices {
    public:
      typedef ExecSpace exec_space;
      typedef ValueType value_type;
      typedef ArrayLayout array_layout;

      typedef BlockTridiagMatrices<exec_space,value_type,array_layout> block_tridiag_matrices_type;

    private:
      ordinal_type _ntridiag, _m, _blocksize, _shmemlvl;

      UnmanagedViewType<typename block_tridiag_matrices_type::value_array_type> _TA, _TB, _TC;
      typedef typename MagnitudeScalarType<value_type>::type magnitude_scalar_type;

    public:
      FactorizeBlockTridiagMatrices() {}

      // A thread maps nonzero blocks
      KOKKOS_INLINE_FUNCTION 
      void operator()(const RangeTag &, const ordinal_type ij) const {
        auto A = Kokkos::subview(_TA, ij, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto B = Kokkos::subview(_TB, ij, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto C = Kokkos::subview(_TC, ij, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());

        auto AA = Kokkos::subview(A, 0,   Kokkos::ALL(), Kokkos::ALL());
        auto BB = Kokkos::subview(B, 0,   Kokkos::ALL(), Kokkos::ALL());
        auto CC = Kokkos::subview(C, 0,   Kokkos::ALL(), Kokkos::ALL());
        auto DD = AA;

        const auto tiny = Kokkos::Details::ArithTraits<magnitude_scalar_type>::epsilon()*100;

        const ordinal_type kend = _m - 1;
        for (ordinal_type k=0;k<kend;++k) {
          AA.assign_data( &A(k  ,0,0) );
          BB.assign_data( &B(k  ,0,0) );
          CC.assign_data( &C(k  ,0,0) );
          DD.assign_data( &A(k+1,0,0) );

          SerialAddRadial::invoke(tiny, AA);
          SerialLU<LU_AlgoTagType>
            ::invoke(AA);
          SerialTrsm<Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::Unit,Trsm_AlgoTagType>
            ::invoke(1.0, AA, BB);
          SerialTrsm<Side::Right,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,Trsm_AlgoTagType>
            ::invoke(1.0, AA, CC);
          SerialGemm<Trans::NoTranspose,Trans::NoTranspose,Gemm_AlgoTagType>
            ::invoke(-1.0, CC, BB, 1.0, DD);
        }
        AA.assign_data( &A(kend,0,0) );
        SerialLU<LU_AlgoTagType>
          ::invoke(AA);
      }

      template<typename MemberType>
      KOKKOS_INLINE_FUNCTION 
      void operator()(const TeamTag &, const MemberType &member) const {
        const int ijbeg = member.league_rank()*VectorLength;
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, VectorLength),
           [&](const int &idx) {
            const int ij = ijbeg + idx;
            if (ij < _ntridiag) {
              auto A = Kokkos::subview(_TA, ij, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
              auto B = Kokkos::subview(_TB, ij, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
              auto C = Kokkos::subview(_TC, ij, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());

              const auto tiny = Kokkos::Details::ArithTraits<magnitude_scalar_type>::epsilon()*100;

              const ordinal_type kend = _m - 1;
              {
                auto AA = Kokkos::subview(A, 0,   Kokkos::ALL(), Kokkos::ALL());
                auto BB = Kokkos::subview(B, 0,   Kokkos::ALL(), Kokkos::ALL());
                auto CC = Kokkos::subview(C, 0,   Kokkos::ALL(), Kokkos::ALL());
                auto DD = AA;

                for (ordinal_type k=0;k<kend;++k) {

                  AA.assign_data( &A(k  ,0,0) );
                  BB.assign_data( &B(k  ,0,0) );
                  CC.assign_data( &C(k  ,0,0) );
                  DD.assign_data( &A(k+1,0,0) );
                  
                  member.team_barrier();
                  TeamAddRadial<MemberType>::invoke(member, tiny, AA);
                  member.team_barrier();
                  TeamLU<MemberType,LU_AlgoTagType>
                    ::invoke(member, AA);
                  member.team_barrier();
                  TeamTrsm<MemberType,Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::Unit,Trsm_AlgoTagType>
                    ::invoke(member, 1.0, AA, BB);
                  TeamTrsm<MemberType,Side::Right,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,Trsm_AlgoTagType>
                    ::invoke(member, 1.0, AA, CC);
                  member.team_barrier();
                  TeamGemm<MemberType,Trans::NoTranspose,Trans::NoTranspose,Gemm_AlgoTagType>
                    ::invoke(member, -1.0, CC, BB, 1.0, DD);
                }
                {
                  member.team_barrier();
                  AA.assign_data( &A(kend,0,0) );
                  TeamLU<MemberType,LU_AlgoTagType>
                    ::invoke(member, AA);
                }
              }
#if 0
              { // 0.028 vs 0.035; without subview it performs 0.028
                const int as0 = A.stride_1(), as1 = A.stride_2();
                const int bs0 = B.stride_1(), bs1 = B.stride_2();
                const int cs0 = C.stride_1(), cs1 = C.stride_2();

                for (ordinal_type k=0;k<kend;++k) {
                  auto AA = &A(k,  0,0);
                  auto BB = &B(k,  0,0);
                  auto CC = &C(k,  0,0);
                  auto DD = &A(k+1,0,0);
                  
                  member.team_barrier();
                  TeamLU_Internal<LU_AlgoTagType>
                    ::invoke(member, _blocksize, _blocksize, AA, as0, as1);
                  member.team_barrier();
                  TeamTrsmInternalLeftLower<Trsm_AlgoTagType>
                    ::invoke(member, true, _blocksize, _blocksize,
                             1.0, AA, as0, as1, BB, bs0, bs1);
                  TeamTrsmInternalLeftLower<Trsm_AlgoTagType>
                    ::invoke(member, false, _blocksize, _blocksize,
                             1.0, AA, as1, as0, CC, cs1, cs0);
                  member.team_barrier();
                  TeamGemmInternal<Gemm_AlgoTagType>::
                    invoke(member, _blocksize, _blocksize, _blocksize, 
                           -1.0, 
                           CC, cs0, cs1, BB, bs0, bs1,
                           1.0,
                           DD, as0, as1);
                }
                {
                  member.team_barrier();
                  auto AA = &A(kend, 0,0);
                  TeamLU_Internal<LU_AlgoTagType>
                    ::invoke(member, _blocksize, _blocksize, AA, as0, as1);
                }
              }
#endif
            }
          });
      }

      template<typename MemberType>
      KOKKOS_INLINE_FUNCTION 
      void operator()(const TeamShmemTag &, const MemberType &member) const {
        typedef Kokkos::View<ValueType***,exec_space> packed_view_type;
        ScratchViewType<packed_view_type> sA(member.team_scratch(_shmemlvl), VectorLength, _blocksize, _blocksize);

        const int ijbeg = member.league_rank()*VectorLength;
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, VectorLength),
           [&](const int &idx) {
            const int ij = ijbeg + idx;
            if (ij < _ntridiag) {
              auto A = Kokkos::subview(_TA, ij, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
              auto B = Kokkos::subview(_TB, ij, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
              auto C = Kokkos::subview(_TC, ij, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());

              auto sAA = Kokkos::subview(sA, idx, Kokkos::ALL(), Kokkos::ALL());
              const ordinal_type kend = _m - 1;
              for (ordinal_type k=0;k<kend;++k) {
                auto AA = Kokkos::subview(A, k,   Kokkos::ALL(), Kokkos::ALL());
                auto BB = Kokkos::subview(B, k,   Kokkos::ALL(), Kokkos::ALL());
                auto CC = Kokkos::subview(C, k,   Kokkos::ALL(), Kokkos::ALL());
                auto DD = Kokkos::subview(A, k+1, Kokkos::ALL(), Kokkos::ALL());

                TeamCopy<MemberType,Trans::NoTranspose>::invoke(member, AA, sAA);
                member.team_barrier();

                TeamLU<MemberType,LU_AlgoTagType>
                  ::invoke(member, sAA);
                member.team_barrier();

                TeamCopy<MemberType,Trans::NoTranspose>::invoke(member, sAA, AA);

                TeamTrsm<MemberType,Side::Left,Uplo::Lower,Trans::NoTranspose,Diag::Unit,Trsm_AlgoTagType>
                  ::invoke(member, 1.0, sAA, BB);
                TeamTrsm<MemberType,Side::Right,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,Trsm_AlgoTagType>
                  ::invoke(member, 1.0, sAA, CC);
                member.team_barrier();

                TeamGemm<MemberType,Trans::NoTranspose,Trans::NoTranspose,Gemm_AlgoTagType>
                  ::invoke(member, -1.0, CC, BB, 1.0, DD);
              }

              {
                member.team_barrier();
                auto AA = Kokkos::subview(A, kend, Kokkos::ALL(), Kokkos::ALL());
                TeamLU<MemberType,LU_AlgoTagType>
                  ::invoke(member, AA);
              }
            }
          });
      }
      
      double FlopCount(const block_tridiag_matrices_type T) {
        const int 
          ntridiag = T.NumTridiagMatrices(),
          m = T.NumRows(),
          blocksize = T.BlockSize();
        
        return ntridiag*( (m-1)*(LU_FlopCount(blocksize, blocksize) + 
                                 Trsm_Lower_FlopCountLower(blocksize, blocksize) +
                                 Trsm_Upper_FlopCountUpper(blocksize, blocksize) +
                                 Gemm_FlopCount(blocksize, blocksize, blocksize)) +
                          LU_FlopCount(blocksize, blocksize) );
      }

      // for batched blas check
      void run(const int op, const block_tridiag_matrices_type T, const bool fake = false) { 
        _ntridiag = T.NumTridiagMatrices();
        _m = T.NumRows(); 
        _blocksize = T.BlockSize();

        _TA = T.A(); 
        _TB = T.B(); 
        _TC = T.C();

        // parallel over the instances of tridiagonal matrices
        if (!fake) {
#if defined(KOKKOS_ENABLE_CUDA) && defined(__KOKKOSBATCHED_TEST_ENABLE_CUDA__)
          typedef FactorizeBlockTridiagMatrices<exec_space,value_type,array_layout,
	    VectorLength,
	    LU_AlgoTagType,
	    Trsm_AlgoTagType,
	    Gemm_AlgoTagType> functor_type;
#endif

          switch (op) {
          case 0: {
            std::cout << "KokkosBatched::RangeTag::" << Gemm_AlgoTagType::name() << "\n";
            const Kokkos::RangePolicy<exec_space,RangeTag> policy(0, _ntridiag);
            Kokkos::parallel_for("KokkosBatched::Test::BlockCrs::FactorizeBlockTridiagMatrices::Op0", policy, *this);
            break;
          }
#if defined(KOKKOS_ENABLE_CUDA) && defined(__KOKKOSBATCHED_TEST_ENABLE_CUDA__)
          case 1: {
            typedef Kokkos::TeamPolicy<exec_space,TeamTag> policy_type;

            int team_size = 0;

            // this is what cuda allows
            const int max_team_size = 
              policy_type(_ntridiag, Kokkos::AUTO, VectorLength).team_size_max(functor_type(), Kokkos::ParallelForTag());

            // this is what algorithm allows
            if (std::is_same<Gemm_AlgoTagType,Algo::Gemm::Blocked>::value) {
              const int mb = Algo::Gemm::Blocked::mb<typename exec_space::memory_space>();
              const int mp = _blocksize%mb, mblk = (_blocksize/mb) + (mp>0);
              // - max parallelism in gemm / 2 (no idea...)
              team_size = std::min(std::max(mblk*mblk/2,1), max_team_size/2);
              std::cout << "KokkosBatched::TeamTag::Blocked::TeamSize:: " << team_size << " " << (max_team_size) << "\n";
            } else {
              // - max parallelism in trsm * scheduling efficiency 2
              // - max cuda team size / scheduling efficiency 2
              team_size = std::min(std::max(_blocksize*2,4), max_team_size/2);
              std::cout << "KokkosBatched::TeamTag::Unblocked::TeamSize:: " << team_size << " " << (max_team_size) << "\n";
            }

            const policy_type policy(_ntridiag, team_size, VectorLength);
            Kokkos::parallel_for("KokkosBatched::Test::BlockCrs::FactorizeBlockTridiagMatrices::Op1", policy, *this);
            break;
          }
          case 2: {
            typedef Kokkos::View<ValueType***,exec_space> packed_view_type;
            typedef Kokkos::TeamPolicy<exec_space,TeamShmemTag> policy_type;
            
            const int per_team_scratch 
              = ScratchViewType<packed_view_type>::shmem_size(VectorLength, _blocksize, _blocksize);
            
            _shmemlvl = ((per_team_scratch/1024) < 48 ? 0 : 1);
            {
              int team_size = 0;
              
              // this is what cuda allows
              const int max_team_size =
                policy_type(_ntridiag, Kokkos::AUTO, VectorLength).set_scratch_size(_shmemlvl, Kokkos::PerTeam(per_team_scratch))
                .team_size_max(functor_type(), Kokkos::ParallelForTag());   
              
              // this is what algorithm allows
              if (std::is_same<Gemm_AlgoTagType,Algo::Gemm::Blocked>::value) {
                const int mb = Algo::Gemm::Blocked::mb<typename exec_space::memory_space>();
                const int mp = _blocksize%mb, mblk = (_blocksize/mb) + (mp>0);
                // - max parallelism in gemm / 2 (no idea...)
                team_size = std::min(std::max(mblk*mblk/2,1), max_team_size/2);
                std::cout << "KokkosBatched::TeamShmemTag::Blocked::TeamSize:: " << team_size << " " << (max_team_size) << "\n";
              } else {
                team_size = std::min(std::max(_blocksize*2,4), max_team_size/2);
                std::cout << "KokkosBatched::TeamShmemTag::Unblocked::TeamSize:: " << team_size << " " << (max_team_size) << "\n";
              }
              
              policy_type policy = policy_type(_ntridiag, team_size, VectorLength).set_scratch_size(_shmemlvl, Kokkos::PerTeam(per_team_scratch));
              Kokkos::parallel_for("KokkosBatched::Test::BlockCrs::FactorizeBlockTridiagMatrices::Op2", policy, *this);
            } 
            break;
          }
#endif
          default: {
            std::cout << "Not supported operation mode: " << op << " \n";
            break;
          }
          }
          
        }
      }

      template<typename AViewType, typename LViewType, typename UViewType>
      void a_subtract_mult_l_and_u(const ordinal_type tl, const ordinal_type il, AViewType A,
                                   const ordinal_type tr, const ordinal_type ir, LViewType L, UViewType U) {
        for (ordinal_type ii=0;ii<_blocksize;++ii)
          for (ordinal_type jj=0;jj<_blocksize;++jj)
            for (ordinal_type kk=0;kk<_blocksize;++kk) {
              const auto l = ( ii == kk ? 1 :
                               ii >  kk ? tdiag_val(L, tr, ir, ii, kk) : 0 );
              const auto u = ( kk <= jj ? tdiag_val(U, tr, ir, kk, jj) : 0 );
              tdiag_val(A, tl, il, ii, jj) -= l*u;
            }
      }

      template<typename AViewType, typename BViewType, typename CViewType>
      void a_subtract_mult_b_and_c(const ordinal_type tl, const ordinal_type il, AViewType A,
                                   const ordinal_type tr, const ordinal_type ir, BViewType B, CViewType C) {
        for (ordinal_type ii=0;ii<_blocksize;++ii)
          for (ordinal_type jj=0;jj<_blocksize;++jj)
            for (ordinal_type kk=0;kk<_blocksize;++kk) 
              tdiag_val(A, tl, il, ii, jj) -= ( tdiag_val(B, tr, ir, ii, kk)*
                                                tdiag_val(C, tr, ir, kk, jj) );
      }

      template<typename AViewType, typename LViewType, typename BViewType>
      void a_subtract_mult_l_and_b(const ordinal_type tl, const ordinal_type il, AViewType A,
                                   const ordinal_type tr, const ordinal_type ir, LViewType L, BViewType B) {
        for (ordinal_type ii=0;ii<_blocksize;++ii)
          for (ordinal_type jj=0;jj<_blocksize;++jj)
            for (ordinal_type kk=0;kk<_blocksize;++kk) {
              const auto l = ( ii == kk ? 1.0 :
                               ii >  kk ? tdiag_val(L, tr, ir, ii, kk) : 0 );
              tdiag_val(A, tl, il, ii, jj) -= l*tdiag_val(B, tr, ir, kk, jj);
            }
      }
      
      template<typename AViewType, typename BViewType, typename UViewType>
      void a_subtract_mult_b_and_u(const ordinal_type tl, const ordinal_type il, AViewType A,
                                   const ordinal_type tr, const ordinal_type ir, BViewType B, UViewType U) {
        for (ordinal_type ii=0;ii<_blocksize;++ii)
          for (ordinal_type jj=0;jj<_blocksize;++jj)
            for (ordinal_type kk=0;kk<_blocksize;++kk) {
              const auto u = ( kk <= jj ? tdiag_val(U, tr, ir, kk, jj) : 0 );
              tdiag_val(A, tl, il, ii, jj) -= tdiag_val(B, tr, ir, ii, kk)*u;
            }
      }

      bool check(const block_tridiag_matrices_type T) {
        // factors
        auto DD = Kokkos::create_mirror_view(_TA);    Kokkos::deep_copy(DD, _TA);
        auto UU = Kokkos::create_mirror_view(_TB);    Kokkos::deep_copy(UU, _TB);
        auto LL = Kokkos::create_mirror_view(_TC);    Kokkos::deep_copy(LL, _TC);

        // input A
        auto A = Kokkos::create_mirror_view(T.A());   Kokkos::deep_copy(A, T.A());
        auto B = Kokkos::create_mirror_view(T.B());   Kokkos::deep_copy(B, T.B());
        auto C = Kokkos::create_mirror_view(T.C());   Kokkos::deep_copy(C, T.C());
        
        // diffs
        Kokkos::View<value_type****,Kokkos::DefaultHostExecutionSpace> 
          AA("AA", _ntridiag, _m,   _blocksize, _blocksize),
          BB("BB", _ntridiag, _m-1, _blocksize, _blocksize),
          CC("CC", _ntridiag, _m-1, _blocksize, _blocksize);
        
        Kokkos::deep_copy(AA, A);
        Kokkos::deep_copy(BB, B);
        Kokkos::deep_copy(CC, C);

        // Check | A - L U | / | A |
        for (ordinal_type t=0;t<_ntridiag;++t) {
          a_subtract_mult_l_and_u(t, 0, AA, 
                                  t, 0, DD, DD);
          for (ordinal_type i=1;i<_m;++i) {
            a_subtract_mult_l_and_u(t, i,   AA, 
                                    t, i,   DD, DD);
            a_subtract_mult_b_and_c(t, i,   AA, 
                                    t, i-1, LL, UU);
            a_subtract_mult_l_and_b(t, i-1, BB, 
                                    t, i-1, DD, UU);
            a_subtract_mult_b_and_u(t, i-1, CC, 
                                    t, i-1, LL, DD);
          }
        }

        double norm = 0, diff = 0;
        for (ordinal_type t=0;t<_ntridiag;++t) {
          for (ordinal_type ii=0;ii<_blocksize;++ii)
            for (ordinal_type jj=0;jj<_blocksize;++jj) {
              norm += std::abs(tdiag_val(A ,t, 0, ii, jj));
              diff += std::abs(tdiag_val(AA,t, 0, ii, jj));
            }
          for (ordinal_type i=1;i<_m;++i) 
            for (ordinal_type ii=0;ii<_blocksize;++ii)
              for (ordinal_type jj=0;jj<_blocksize;++jj) {
                norm += std::abs(tdiag_val(A ,t, i,   ii, jj));
                diff += std::abs(tdiag_val(AA,t, i,   ii, jj));
                norm += std::abs(tdiag_val(B ,t, i-1, ii, jj));
                diff += std::abs(tdiag_val(BB,t, i-1, ii, jj));
                norm += std::abs(tdiag_val(C ,t, i-1, ii, jj));
                diff += std::abs(tdiag_val(CC,t, i-1, ii, jj));
              }
        }
        //std::cout << "tridiag factor check  norm = " << norm << "  diff = " << diff << std::endl;
        const bool r_val = diff/norm < 1e2*std::numeric_limits<scalar_type>::epsilon();
        return r_val;
      }
    };

    template <typename ExecSpace, typename ValueType, typename ArrayLayout,
              int VectorLength,
              typename Trsv_AlgoTagType,
              typename Gemv_AlgoTagType>
    class SolveBlockTridiagMatrices {
    public:
      typedef ExecSpace exec_space;
      typedef ValueType value_type;
      typedef ArrayLayout array_layout;

      typedef BlockTridiagMatrices<exec_space,value_type,array_layout> block_tridiag_matrices_type;
      typedef PartitionedBlockMultiVector<exec_space,value_type,array_layout> partitioned_block_multi_vector_type;
      
    private:
      ordinal_type _ntridiag, _m, _blocksize, _nvectors, _shmemlvl;

      ConstUnmanagedViewType<typename block_tridiag_matrices_type::value_array_type> _TA, _TB, _TC;
      ConstUnmanagedViewType<typename partitioned_block_multi_vector_type::value_array_type> _b;
      /**/ UnmanagedViewType<typename partitioned_block_multi_vector_type::value_array_type> _x;

    public:
      SolveBlockTridiagMatrices() {}

      KOKKOS_INLINE_FUNCTION 
      void operator()(const RangeTag &, const ordinal_type ij) const {
        auto A = Kokkos::subview(_TA, ij, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto B = Kokkos::subview(_TB, ij, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
        auto C = Kokkos::subview(_TC, ij, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());

        // subview patterns
        auto A_0_all_all = Kokkos::subview(A, 0, Kokkos::ALL(), Kokkos::ALL());
        auto B_0_all_all = Kokkos::subview(B, 0, Kokkos::ALL(), Kokkos::ALL());
        auto C_0_all_all = Kokkos::subview(C, 0, Kokkos::ALL(), Kokkos::ALL());
        
        auto mx_0_0_all_all = Kokkos::subview(_x, 0, 0, Kokkos::ALL(), Kokkos::ALL());
        auto mb_0_0_all_all = Kokkos::subview(_b, 0, 0, Kokkos::ALL(), Kokkos::ALL());

        auto x_0_all = Kokkos::subview(mx_0_0_all_all, 0, Kokkos::ALL());
        auto b_0_all = Kokkos::subview(mb_0_0_all_all, 0, Kokkos::ALL());

        ///
        /// loop over multivectors
        ///

        auto &x = mx_0_0_all_all;
        auto &b = mb_0_0_all_all;

        auto &xt = x_0_all;
        auto  xb = x_0_all;

        auto &bt = b_0_all;
        auto  bb = b_0_all;

        for (int jvec=0;jvec<_nvectors;++jvec) {
          x.assign_data( &_x(ij,jvec,0,0) );
          b.assign_data( &_b(ij,jvec,0,0) );
          
          ///
          /// forward substitution
          ///
          {
            auto &LT = A_0_all_all; 
            auto &LB = C_0_all_all; 
            
            const bool is_same_x_and_b = (x.data() == b.data());
            {
              if (!is_same_x_and_b) {
                xt.assign_data( &x(0,0) );
                bt.assign_data( &b(0,0) );
                SerialCopy<Trans::NoTranspose>::invoke(bt, xt);
              }
            }
            const ordinal_type kend = _m - 1;
            for (ordinal_type k=0;k<kend;++k) {
              LT.assign_data( &A(k,0,0) );
              LB.assign_data( &C(k,0,0) );
              
              xt.assign_data( &x(k  ,0) );
              xb.assign_data( &x(k+1,0) );

              if (!is_same_x_and_b) {
                bb.assign_data( &b(k+1,0) );
                SerialCopy<Trans::NoTranspose>::invoke(bb, xb);
              }
              
              SerialTrsv<Uplo::Lower,Trans::NoTranspose,Diag::Unit,Trsv_AlgoTagType>
                ::invoke(1.0, LT, xt);
              SerialGemv<Trans::NoTranspose,Gemv_AlgoTagType>
                ::invoke(-1.0, LB, xt, 1.0, xb);
            }
            
            LT.assign_data( &A(kend,0,0) );
            xt.assign_data( &x(kend,0) );
            SerialTrsv<Uplo::Lower,Trans::NoTranspose,Diag::Unit,Trsv_AlgoTagType>
              ::invoke(1.0, LT, xt);
          }

          ///
          /// backward substitution
          ///
          {
            auto &UT = B_0_all_all;
            auto &UB = A_0_all_all;

            const ordinal_type kbegin = _m - 1;
            for (ordinal_type k=kbegin;k>0;--k) {
              UT.assign_data( &B(k-1,0,0) );
              UB.assign_data( &A(k  ,0,0) );
              
              xt.assign_data( &x(k-1,0) );
              xb.assign_data( &x(k  ,0) );

              SerialTrsv<Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,Trsv_AlgoTagType>
                ::invoke(1.0, UB, xb);
              SerialGemv<Trans::NoTranspose,Gemv_AlgoTagType>
                ::invoke(-1.0, UT, xb, 1.0, xt);
            }
            UT.assign_data( &A(0,0,0) );
            xt.assign_data( &x(0,0) );
            SerialTrsv<Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,Trsv_AlgoTagType>
              ::invoke(1.0, UT, xt);
          }
        }
      }

      template<typename MemberType>
      KOKKOS_INLINE_FUNCTION 
      void operator()(const TeamTag &, const MemberType &member) const {
        const int ijbeg = member.league_rank()*VectorLength;
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, VectorLength),
           [&](const int &idx) {
            const int ij = ijbeg + idx;
            if (ij < _ntridiag) {
              auto A = Kokkos::subview(_TA, ij, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
              auto B = Kokkos::subview(_TB, ij, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
              auto C = Kokkos::subview(_TC, ij, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
              
              ///
              /// loop over multivectors
              ///
              for (int jvec=0;jvec<_nvectors;++jvec) {
                auto x = Kokkos::subview(_x, ij, jvec, Kokkos::ALL(), Kokkos::ALL());
                auto b = Kokkos::subview(_b, ij, jvec, Kokkos::ALL(), Kokkos::ALL());
                
                ///
                /// forward substitution
                ///
                {
                  const bool is_same_x_and_b = (x.data() == b.data());
                  {
                    if (!is_same_x_and_b) {
                      auto x0 = Kokkos::subview(x, 0, Kokkos::ALL());          
                      auto b0 = Kokkos::subview(b, 0, Kokkos::ALL());
                      TeamCopy<MemberType,Trans::NoTranspose>::invoke(member, b0, x0);
                      member.team_barrier();
                    }
                  }
                  const ordinal_type kend = _m - 1;
                  for (ordinal_type k=0;k<kend;++k) {
                    auto LT = Kokkos::subview(A, k,   Kokkos::ALL(), Kokkos::ALL());
                    auto LB = Kokkos::subview(C, k,   Kokkos::ALL(), Kokkos::ALL());

                    auto xt = Kokkos::subview(x, k,   Kokkos::ALL());
                    auto xb = Kokkos::subview(x, k+1, Kokkos::ALL());

                    if (!is_same_x_and_b) {
                      auto bb = Kokkos::subview(b, k+1, Kokkos::ALL());
                      TeamCopy<MemberType,Trans::NoTranspose>::invoke(member, bb, xb);
                    }

                    member.team_barrier();
                    TeamTrsv<MemberType,Uplo::Lower,Trans::NoTranspose,Diag::Unit,Trsv_AlgoTagType>
                      ::invoke(member, 1.0, LT, xt);

                    member.team_barrier();
                    TeamGemv<MemberType,Trans::NoTranspose,Gemv_AlgoTagType>
                      ::invoke(member, -1.0, LB, xt, 1.0, xb);
                  }
                  {
                    auto LL = Kokkos::subview(A, kend, Kokkos::ALL(), Kokkos::ALL());
                    auto xx = Kokkos::subview(x, kend, Kokkos::ALL());
                    member.team_barrier();
                    TeamTrsv<MemberType,Uplo::Lower,Trans::NoTranspose,Diag::Unit,Trsv_AlgoTagType>
                      ::invoke(member, 1.0, LL, xx);
                  }
                }

                ///
                /// backward substitution
                ///
                {
                  const ordinal_type kbegin = _m - 1;
                  for (ordinal_type k=kbegin;k>0;--k) {
                    auto UT = Kokkos::subview(B, k-1, Kokkos::ALL(), Kokkos::ALL());
                    auto UB = Kokkos::subview(A, k,   Kokkos::ALL(), Kokkos::ALL());
            
                    auto xt = Kokkos::subview(x, k-1, Kokkos::ALL());
                    auto xb = Kokkos::subview(x, k,   Kokkos::ALL());

                    member.team_barrier();
                    TeamTrsv<MemberType,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,Trsv_AlgoTagType>
                      ::invoke(member, 1.0, UB, xb);

                    member.team_barrier();
                    TeamGemv<MemberType,Trans::NoTranspose,Gemv_AlgoTagType>
                      ::invoke(member, -1.0, UT, xb, 1.0, xt);
                  }
                  {
                    auto UU = Kokkos::subview(A, 0, Kokkos::ALL(), Kokkos::ALL());
                    auto xx = Kokkos::subview(x, 0, Kokkos::ALL());

                    member.team_barrier();
                    TeamTrsv<MemberType,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,Trsv_AlgoTagType>
                      ::invoke(member, 1.0, UU, xx);
                  }
                }
              }
            }
          });
      }
      
      template<typename MemberType>
      KOKKOS_INLINE_FUNCTION 
      void operator()(const TeamShmemTag &, const MemberType &member) const {
        typedef Kokkos::View<ValueType***,exec_space> packed_view_type;        
        ScratchViewType<packed_view_type> s(member.team_scratch(_shmemlvl), VectorLength, _m, _blocksize);

        const int ijbeg = member.league_rank()*VectorLength;
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(member, VectorLength),
           [&](const int &idx) {
            const int ij = ijbeg + idx;
            if (ij < _ntridiag) {
              auto A = Kokkos::subview(_TA, ij, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
              auto B = Kokkos::subview(_TB, ij, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());
              auto C = Kokkos::subview(_TC, ij, Kokkos::ALL(), Kokkos::ALL(), Kokkos::ALL());

              auto sx = Kokkos::subview(s, idx, Kokkos::ALL(), Kokkos::ALL());              

              ///
              /// loop over multivectors
              ///
              for (int jvec=0;jvec<_nvectors;++jvec) {
                auto x = Kokkos::subview(_x, ij, jvec, Kokkos::ALL(), Kokkos::ALL());
                auto b = Kokkos::subview(_b, ij, jvec, Kokkos::ALL(), Kokkos::ALL());

                // copy the entire vector into shared memory (if necessary it needs chunking)
                TeamCopy<MemberType,Trans::NoTranspose>::invoke(member, b, sx);
                member.team_barrier();
                
                ///
                /// forward substitution
                ///
                {
                  const ordinal_type kend = _m - 1;
                  for (ordinal_type k=0;k<kend;++k) {
                    auto LT = Kokkos::subview(A, k,   Kokkos::ALL(), Kokkos::ALL());
                    auto LB = Kokkos::subview(C, k,   Kokkos::ALL(), Kokkos::ALL());

                    auto xt = Kokkos::subview(sx, k,   Kokkos::ALL());
                    auto xb = Kokkos::subview(sx, k+1, Kokkos::ALL());

                    member.team_barrier();
                    TeamTrsv<MemberType,Uplo::Lower,Trans::NoTranspose,Diag::Unit,Trsv_AlgoTagType>
                      ::invoke(member, 1.0, LT, xt);

                    member.team_barrier();
                    TeamGemv<MemberType,Trans::NoTranspose,Gemv_AlgoTagType>
                      ::invoke(member, -1.0, LB, xt, 1.0, xb);
                  }
                  {
                    auto LL = Kokkos::subview(A, kend, Kokkos::ALL(), Kokkos::ALL());
                    auto xx = Kokkos::subview(sx, kend, Kokkos::ALL());
                    member.team_barrier();
                    TeamTrsv<MemberType,Uplo::Lower,Trans::NoTranspose,Diag::Unit,Trsv_AlgoTagType>
                      ::invoke(member, 1.0, LL, xx);
                  }
                }

                ///
                /// backward substitution
                ///
                {
                  const ordinal_type kbegin = _m - 1;
                  for (ordinal_type k=kbegin;k>0;--k) {
                    auto UT = Kokkos::subview(B, k-1, Kokkos::ALL(), Kokkos::ALL());
                    auto UB = Kokkos::subview(A, k,   Kokkos::ALL(), Kokkos::ALL());
            
                    auto xt = Kokkos::subview(sx, k-1, Kokkos::ALL());
                    auto xb = Kokkos::subview(sx, k,   Kokkos::ALL());

                    member.team_barrier();
                    TeamTrsv<MemberType,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,Trsv_AlgoTagType>
                      ::invoke(member, 1.0, UB, xb);

                    member.team_barrier();
                    TeamGemv<MemberType,Trans::NoTranspose,Gemv_AlgoTagType>
                      ::invoke(member, -1.0, UT, xb, 1.0, xt);
                  }
                  {
                    auto UU = Kokkos::subview(A, 0, Kokkos::ALL(), Kokkos::ALL());
                    auto xx = Kokkos::subview(sx, 0, Kokkos::ALL());

                    member.team_barrier();
                    TeamTrsv<MemberType,Uplo::Upper,Trans::NoTranspose,Diag::NonUnit,Trsv_AlgoTagType>
                      ::invoke(member, 1.0, UU, xx);
                  }
                }
                // copy the entire vector into shared memory (if necessary it needs chunking)
                TeamCopy<MemberType,Trans::NoTranspose>::invoke(member, sx, x);
                member.team_barrier();
              }
            }
          });


      }

      void run(const int op, 
               const block_tridiag_matrices_type T,
               const partitioned_block_multi_vector_type x,
               const partitioned_block_multi_vector_type b) {
        assert(T.NumTridiagMatrices() == x.NumPartitions());
        assert(T.NumRows() == x.NumRows());
        assert(T.BlockSize() == x.BlockSize());

        _ntridiag = T.NumTridiagMatrices();
        _m = T.NumRows(); 
        _blocksize = T.BlockSize();
        _nvectors = x.NumVectors();

        _TA = T.A(); 
        _TB = T.B(); 
        _TC = T.C();
        
        _x = x.Values();
        _b = b.Values();

        {
#if defined(KOKKOS_ENABLE_CUDA) && defined(__KOKKOSBATCHED_TEST_ENABLE_CUDA__)
          typedef SolveBlockTridiagMatrices<exec_space,value_type,array_layout,
	    VectorLength,
	    Trsv_AlgoTagType,
	    Gemv_AlgoTagType> functor_type;
#endif

          switch (op) {
          case 0: {
            const Kokkos::RangePolicy<exec_space,RangeTag> policy(0, _ntridiag);
            Kokkos::parallel_for("KokkosBatched::Test::BlockCrs::SolveBlockTridiagMatrices::Op0", policy, *this);
            break;
          }
#if defined(KOKKOS_ENABLE_CUDA) && defined(__KOKKOSBATCHED_TEST_ENABLE_CUDA__)
          case 1: {
            typedef Kokkos::TeamPolicy<exec_space,TeamTag> policy_type;
            
            int team_size = 0;
            
            // this is what cuda allows
            const int max_team_size = 
              policy_type(_ntridiag, Kokkos::AUTO, VectorLength).team_size_max(functor_type(), Kokkos::ParallelForTag());

            // this is what algorithm allows
            if (std::is_same<Gemv_AlgoTagType,Algo::Gemv::Blocked>::value) {
              const int mb = Algo::Gemv::Blocked::mb<typename exec_space::memory_space>();
              const int mp = _blocksize%mb, mblk = (_blocksize/mb) + (mp>0);
              team_size = std::min(std::max(mblk/2,1), int(max_team_size/2));
            } else {
              // in solve phase, max peak parallelism is same as blocksize (one iteration)
              // better to give blocksize/2 
              team_size = std::min(std::max(_blocksize/2,4), int(max_team_size/2));
            }
            
            const policy_type policy(_ntridiag, team_size, VectorLength);
            Kokkos::parallel_for("KokkosBatched::Test::BlockCrs::SolveBlockTridiagMatrices::Op1", policy, *this);
            break;
          }
          case 2: {
            typedef Kokkos::View<ValueType***,exec_space> packed_view_type;
            typedef Kokkos::TeamPolicy<exec_space,TeamShmemTag> policy_type;

            const int per_team_scratch 
              = ScratchViewType<packed_view_type>::shmem_size(VectorLength, _m, _blocksize);

            _shmemlvl = ((per_team_scratch/1024) < 48 ? 0 : 1);            
            {
              int team_size = 0;
              
              // this is what cuda allows
              const int max_team_size =
                policy_type(_ntridiag, Kokkos::AUTO, VectorLength).set_scratch_size(_shmemlvl, Kokkos::PerTeam(per_team_scratch))
                .team_size_max(functor_type(), Kokkos::ParallelForTag());   
              
              // this is what algorithm allows
              if (std::is_same<Gemv_AlgoTagType,Algo::Gemv::Blocked>::value) {
                const int mb = Algo::Gemv::Blocked::mb<typename exec_space::memory_space>();
                const int mp = _blocksize%mb, mblk = (_blocksize/mb) + (mp>0);
                team_size = std::min(std::max(mblk/2,1), int(max_team_size/2));
              } else {
                team_size = std::min(std::max(_blocksize/2,4), int(max_team_size/2));
              }
              
              policy_type policy = policy_type(_ntridiag, team_size, VectorLength).set_scratch_size(_shmemlvl, Kokkos::PerTeam(per_team_scratch));;
              Kokkos::parallel_for("KokkosBatched::Test::BlockCrs::SolveBlockTridiagMatrices::Op2", policy, *this);
            }
            break;
          }
#endif
          default: {
            std::cout << "Not supported operation mode: " << op << " \n";
            break;
          }
          }
        }
      }

      template<typename RViewType, typename AViewType, typename XViewType>
      void r_subtract_mult_a_and_x(const ordinal_type tr, const ordinal_type ir, RViewType R,
                                   const ordinal_type ta, const ordinal_type ia, AViewType A,
                                   const ordinal_type tx, const ordinal_type ix, XViewType X) {
        for (ordinal_type kk=0;kk<_nvectors;++kk) 
          for (ordinal_type ii=0;ii<_blocksize;++ii)
            for (ordinal_type jj=0;jj<_blocksize;++jj)
              tdiag_val(R, tr, kk, ir, ii) -= tdiag_val(A, ta, ia, ii, jj) * tdiag_val(X, tx, kk, ix, jj);
      }

      bool check(const block_tridiag_matrices_type T, 
                 const partitioned_block_multi_vector_type b) {
        // input A
        auto AA = Kokkos::create_mirror_view(T.A());      Kokkos::deep_copy(AA, T.A());
        auto BB = Kokkos::create_mirror_view(T.B());      Kokkos::deep_copy(BB, T.B());
        auto CC = Kokkos::create_mirror_view(T.C());      Kokkos::deep_copy(CC, T.C());

        auto bb = Kokkos::create_mirror_view(b.Values()); Kokkos::deep_copy(bb, b.Values());
        auto xx = Kokkos::create_mirror_view(_x);         Kokkos::deep_copy(xx, _x);
        
        // diffs
        Kokkos::View<value_type****,Kokkos::DefaultHostExecutionSpace> 
          rr("rr", bb.extent(0), bb.extent(1), bb.extent(2), bb.extent(3)); 
        
        Kokkos::deep_copy(rr, bb);

        // Check | Ax - b | / | b |
        for (ordinal_type t=0;t<_ntridiag;++t) {
          r_subtract_mult_a_and_x(t, 0, rr, 
                                  t, 0, AA, 
                                  t, 0, xx);
          r_subtract_mult_a_and_x(t, 0, rr, 
                                  t, 0, BB, 
                                  t, 1, xx);
          
          for (ordinal_type i=1;i<(_m-1);++i) {
            r_subtract_mult_a_and_x(t, i,   rr, 
                                    t, i-1, CC, 
                                    t, i-1, xx);
            r_subtract_mult_a_and_x(t, i,   rr, 
                                    t, i,   AA, 
                                    t, i,   xx);
            r_subtract_mult_a_and_x(t, i,   rr, 
                                    t, i,   BB, 
                                    t, i+1, xx);
          }
          r_subtract_mult_a_and_x(t, _m-1, rr, 
                                  t, _m-2, CC, 
                                  t, _m-2, xx);
          r_subtract_mult_a_and_x(t, _m-1, rr, 
                                  t, _m-1, AA, 
                                  t, _m-1, xx);
        }

        double norm = 0, diff = 0;
        for (ordinal_type t=0;t<_ntridiag;++t) 
          for (ordinal_type jvec=0;jvec<_nvectors;++jvec) 
            for (ordinal_type i=0;i<_m;++i) 
              for (ordinal_type ii=0;ii<_blocksize;++ii) {
                norm += std::abs(tdiag_val(bb, t, jvec, i, ii));
                diff += std::abs(tdiag_val(rr, t, jvec, i, ii));
              }

        //std::cout << "tridiag solve check  norm = " << norm << "  diff = " << diff << std::endl;
        const bool r_val = diff/norm < 1e2*std::numeric_limits<scalar_type>::epsilon();
        return r_val;
      }
    };

    // unit tests
    template<typename DeviceSpace, typename ValueType, int VectorLength, int Oper>
    void run(const ordinal_type ni, const ordinal_type nj, const ordinal_type nk, 
             const ordinal_type blocksize, 
             const ordinal_type nrhs,
             const bool test_tpl = false) {
      typedef typename DeviceSpace::array_layout DeviceArrayLayout;
      typedef Kokkos::DefaultHostExecutionSpace HostSpace;

      bool success = true;
      StructuredBlock mesh(ni, nj, nk);

      // Test StructuredBlock.
      for (ordinal_type c=0;c<mesh.size();++c) {
        ordinal_type i, j, k;
        mesh.id2ijk(c, i, j, k);
        TEST_ASSERT(i >= 0 && i < mesh.ni, success);
        TEST_ASSERT(j >= 0 && j < mesh.nj, success);
        TEST_ASSERT(k >= 0 && k < mesh.nk, success);
        TEST_ASSERT(mesh.ijk2id(i, j, k) == c, success);
      }

      // Graph construction
      CrsGraph<HostSpace,DeviceArrayLayout> graph_host 
        = create_graph_host_for_structured_block<DeviceArrayLayout>(mesh, StencilShape::cross);

      // Crs matrix and multi vector construction
      BlockCrsMatrix<HostSpace,DeviceArrayLayout> A_host(graph_host, blocksize);
      fill_block_crs_matrix_host(A_host);      

      // Device mirroring
      auto A_device = create_mirror<DeviceSpace>(A_host);
      deep_copy(A_device, A_host);

      // Test Matrix Vector product
      {
        const ordinal_type m = graph_host.NumRows();
        
        BlockMultiVector<HostSpace,DeviceArrayLayout> x_host(nrhs, m, blocksize);
        fill_block_multi_vector_host(x_host);
        
        auto x_device = create_mirror<DeviceSpace>(x_host);
        deep_copy(x_device, x_host);
        
        BlockMultiVector<DeviceSpace,DeviceArrayLayout> 
          y1_device(nrhs, m, blocksize), 
          y2_device(nrhs, m, blocksize);

        {
          BlockCrsMatrixVectorProductByRow<DeviceSpace,DeviceArrayLayout> matvec;
          matvec.run(A_device, x_device, y1_device);
        }
        {
          BlockCrsMatrixVectorProductByBlockRow<DeviceSpace,DeviceArrayLayout> matvec;
          matvec.run(A_device, x_device, y2_device);
        }

        const double rdiff = compute_relative_diff(y1_device.Values(), y2_device.Values());
        TEST_ASSERT(rdiff <= 1e2*std::numeric_limits<scalar_type>::epsilon(), success);
      }

      // Test Block TriDiag Extraction
      BlockTridiagMatrices<DeviceSpace,ValueType,DeviceArrayLayout> T_device
        = create_block_tridiag_matrices
        <DeviceSpace,ValueType,DeviceArrayLayout>(mesh.ni*mesh.nj, 
                                                  mesh.nk,
                                                  blocksize);
      {
        ExtractBlockTridiagMatrices<DeviceSpace,ValueType,DeviceArrayLayout> extblk(mesh);
        extblk.run(A_device, T_device);
        TEST_ASSERT(extblk.check(), success);
      }

      BlockTridiagMatrices<DeviceSpace,ValueType,DeviceArrayLayout> T_org_device
        = create_block_tridiag_matrices
        <DeviceSpace,ValueType,DeviceArrayLayout>(mesh.ni*mesh.nj, 
                                                  mesh.nk,
                                                  blocksize);
      
      deep_copy(T_org_device, T_device);
      
      // Test Block TriDiag Factorization
      if (test_tpl) {
#if                                                     \
  defined(__KOKKOSBATCHED_INTEL_MKL__) &&               \
  defined(__KOKKOSBATCHED_INTEL_MKL_BATCHED__) &&       \
  defined(__KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__)
        FactorizeBlockTridiagMatrices<DeviceSpace,
	  ValueType,
	  DeviceArrayLayout,
	  VectorLength,
	  Algo::LU::CompactMKL,
	  Algo::Trsm::CompactMKL,
	  Algo::Gemm::CompactMKL> factorblk;
        factorblk.run(0, T_device); // range policy only now
        TEST_ASSERT(factorblk.check(T_org_device), success);
#elif defined(__KOKKOSBATCHED_NVIDIA_CUBLAS__)
        std::cout << "CUBLAS compact version does not exist\n"; 
#else
        std::cout << "TPLs (CompactMKL or CUBLAS) are not found\n"; 
#endif 
      } else {
        FactorizeBlockTridiagMatrices<DeviceSpace,
	  ValueType, 
	  DeviceArrayLayout,
	  VectorLength,
	  AlgoLU,
	  AlgoTrsm,
	  AlgoGemm> factorblk;
        factorblk.run(Oper, T_device);
        TEST_ASSERT(factorblk.check(T_org_device), success);
      }

      // Test Block TriDiag Solve
      {
        PartitionedBlockMultiVector<HostSpace,ValueType,DeviceArrayLayout> b_host
          = create_partitioned_block_multi_vector
          <HostSpace,ValueType,DeviceArrayLayout>(mesh.ni*mesh.nj, 
                                                  nrhs,
                                                  mesh.nk, 
                                                  blocksize);
        fill_partitioned_block_multi_vector_host(b_host, mesh.ni*mesh.nj);

        auto b_device = create_mirror<DeviceSpace>(b_host);
        deep_copy(b_device, b_host);

        PartitionedBlockMultiVector<DeviceSpace,ValueType,DeviceArrayLayout> x_device
          = create_partitioned_block_multi_vector
          <DeviceSpace,ValueType,DeviceArrayLayout>(mesh.ni*mesh.nj, 
                                                    nrhs,
                                                    mesh.nk, 
                                                    blocksize);
        if (test_tpl) {
#if                                                     \
  defined(__KOKKOSBATCHED_INTEL_MKL__) &&               \
  defined(__KOKKOSBATCHED_INTEL_MKL_BATCHED__) &&       \
  defined(__KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__)
          SolveBlockTridiagMatrices<DeviceSpace,
	    ValueType,
	    DeviceArrayLayout,
	    VectorLength,
	    Algo::Trsv::CompactMKL,
	    Algo::Gemv::CompactMKL> solveblk;

          solveblk.run(0, T_device, x_device, b_device);
          TEST_ASSERT(solveblk.check(T_org_device, b_device), success);
#elif defined(__KOKKOSBATCHED_NVIDIA_CUBLAS__)
	  std::cout << "CUBLAS compact version does not exist\n"; 
#else
	  std::cout << "TPLs (CompactMKL or CUBLAS) are not found\n"; 
#endif 
        } else {
          SolveBlockTridiagMatrices<DeviceSpace,
	    ValueType, 
	    DeviceArrayLayout,
	    VectorLength,
	    AlgoTrsv,
	    AlgoGemv> solveblk;
          
          solveblk.run(Oper, T_device, x_device, b_device);
          TEST_ASSERT(solveblk.check(T_org_device, b_device), success);
        }
      }

      if (!success)
        std::cout << "Unit Tests:: Failed:: "
                  << " ni = " << ni << " nj = " << nj << " nk = " << nk 
                  << " blocksize = " << blocksize << " nrhs = " << nrhs << " \n"; 
    }
    
    // performance tests
    template<typename DeviceSpace, typename ValueType, int VectorLength>
    int run(const Input<DeviceSpace> &input, const bool test_tpl = false) { 
      typedef typename DeviceSpace::array_layout DeviceArrayLayout;
      typedef Kokkos::DefaultHostExecutionSpace HostSpace;

      const ordinal_type niter = 50;      
      int dontopt = 0;
      bool success = true;
      
      /// 
      /// construct a discrete system of equations
      ///
      const ordinal_type 
        ni = input.ni, 
        nj = input.nj, 
        nk = input.nk,
        blocksize = input.bs, 
        nrhs = input.nrhs,
        opf =input.opf,
        ops =input.ops;

      StructuredBlock mesh(ni, nj, nk);

      // something is not copyable ... don't know why yet...
      BlockCrsMatrix<DeviceSpace,DeviceArrayLayout> A_device;
      //double t_fill_block_crs_matrix = 0.0, t_fill_graph = 0.0;
      {
        const StencilShape::Enum stencil_shape = input.stencil_shape;
        CrsGraph<HostSpace,DeviceArrayLayout> graph_host;
        {
          Timer timer("Fill Graph _______________");
          timer.reset();
          graph_host = create_graph_host_for_structured_block<DeviceArrayLayout>(mesh, stencil_shape);
          /* t_fill_graph = */ timer.seconds();
        }
        BlockCrsMatrix<HostSpace,DeviceArrayLayout> A_host(graph_host, blocksize);
        {
          Timer timer("Fill Block CRS Matrix_______________");
          timer.reset();
          fill_block_crs_matrix_host(A_host);       
          /* t_fill_block_crs_matrix = */timer.seconds();       
        }
        A_device = create_mirror<DeviceSpace>(A_host);
        deep_copy(A_device, A_host);
      }

      // memory size
      const double memsize_A = A_device.Values().extent(0)*blocksize*blocksize*8;
      
      ///
      /// matrix vector multiplication test
      ///
      double t_matvec = 0.0;
      //double t_fill_block_multi_vector = 0.0;
      {
        const ordinal_type m = mesh.size();

        BlockMultiVector<HostSpace,DeviceArrayLayout> x_host(nrhs, m, blocksize);
        {
          Timer timer("Fill Block Multi Vector______________");
          timer.reset();
          fill_block_multi_vector_host(x_host);
          /* t_fill_block_multi_vector = */ timer.seconds();
        }
        auto x_device = create_mirror<DeviceSpace>(x_host);
        deep_copy(x_device, x_host);
        
        BlockMultiVector<DeviceSpace,DeviceArrayLayout> y_device(nrhs, m, blocksize);
        {
          //BlockCrsMatrixVectorProductByRow<DeviceSpace> matvec;
          BlockCrsMatrixVectorProductByBlockRow<DeviceSpace,DeviceArrayLayout> matvec;
          {
            Timer timer("50 BlockCrsMatrixVectorProduct");
            timer.reset();
            for (ordinal_type i=0;i<niter;++i) {
              matvec.run(A_device, x_device, y_device);
              dontopt += i;
            }
            t_matvec = timer.seconds();
          }
        }
      }

      ///
      /// block tridiag extraction test
      ///
      const double memsize_T = ni*nj*(3*(nk-1)*blocksize*blocksize + blocksize*blocksize)*8;

      double t_extract = 0.0;
      BlockTridiagMatrices<DeviceSpace,ValueType,DeviceArrayLayout> T_device
        = create_block_tridiag_matrices<DeviceSpace,ValueType,DeviceArrayLayout>(ni*nj, nk, blocksize);
      {
        ExtractBlockTridiagMatrices<DeviceSpace,ValueType,DeviceArrayLayout> extblk(mesh);
        {
          Timer timer("ExtractBlockTridiagMatrices");
          timer.reset();
          extblk.run(A_device, T_device);
          t_extract = timer.seconds();
        }
        if (input.check) TEST_ASSERT(extblk.check(), success);
      }

      // keep original matrix for check
      BlockTridiagMatrices<DeviceSpace,ValueType,DeviceArrayLayout> T_org_device
        = create_block_tridiag_matrices<DeviceSpace,ValueType,DeviceArrayLayout>(ni*nj, nk, blocksize);
      
      deep_copy(T_org_device, T_device);
      
      ///
      /// block tridiag factorization test
      ///
      double t_factorize = 0.0, f_factorize = 0.0;
      if (test_tpl) {
#if                                                     \
  defined(__KOKKOSBATCHED_INTEL_MKL__) &&               \
  defined(__KOKKOSBATCHED_INTEL_MKL_BATCHED__) &&       \
  defined(__KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__)
        FactorizeBlockTridiagMatrices<DeviceSpace,
	  ValueType,
	  DeviceArrayLayout,
	  VectorLength,
	  Algo::LU::CompactMKL,
	  Algo::Trsm::CompactMKL,
	  Algo::Gemm::CompactMKL> factorblk;
        f_factorize = factorblk.FlopCount(T_device)*(sizeof(ValueType)/sizeof(double));
        {
          Timer timer("FactorizeBlockTridiagMatrices");
          timer.reset();
          Kokkos::fence();
          factorblk.run(0, T_device);
          Kokkos::fence();
          t_factorize = timer.seconds();
        }
        TEST_ASSERT(factorblk.check(T_org_device), success);
#elif defined(__KOKKOSBATCHED_NVIDIA_CUBLAS__)
        std::cout << "CUBLAS compact version does not exist\n"; 
#else
        std::cout << "TPLs (CompactMKL or CUBLAS) are not found\n"; 
#endif 
      } else {
        FactorizeBlockTridiagMatrices<DeviceSpace,
	  ValueType, 
	  DeviceArrayLayout,
	  VectorLength,
	  AlgoLU,
	  AlgoTrsm,
	  AlgoGemm> factorblk;
        
        f_factorize = factorblk.FlopCount(T_device)*(sizeof(ValueType)/sizeof(double));
        {
          Timer timer("FactorizeBlockTridiagMatrices");
          timer.reset();
          Kokkos::fence();
          factorblk.run(opf, T_device);
          Kokkos::fence();
          t_factorize = timer.seconds();
        }
        if (input.check) TEST_ASSERT(factorblk.check(T_org_device), success);
      }
      
      ///
      /// block tridiag solve test
      ///
      double t_solve = 0.0;
      {
        PartitionedBlockMultiVector<HostSpace,ValueType,DeviceArrayLayout> b_host
          = create_partitioned_block_multi_vector
          <HostSpace,ValueType,DeviceArrayLayout>(ni*nj, 
                                                  nrhs,
                                                  nk, 
                                                  blocksize);
        fill_partitioned_block_multi_vector_host(b_host, ni*nj);

        auto b_device = create_mirror<DeviceSpace>(b_host);
        deep_copy(b_device, b_host);
        
        PartitionedBlockMultiVector<DeviceSpace,ValueType,DeviceArrayLayout> x_device
          = create_partitioned_block_multi_vector
          <DeviceSpace,ValueType,DeviceArrayLayout>(ni*nj, 
                                                    nrhs,
                                                    nk, 
                                                    blocksize);
        if (test_tpl) {
#if                                                     \
  defined(__KOKKOSBATCHED_INTEL_MKL__) &&               \
  defined(__KOKKOSBATCHED_INTEL_MKL_BATCHED__) &&       \
  defined(__KOKKOSBATCHED_INTEL_MKL_COMPACT_BATCHED__)
          SolveBlockTridiagMatrices<DeviceSpace,
	    ValueType, 
	    DeviceArrayLayout,
	    VectorLength,
	    Algo::Trsv::CompactMKL,
	    Algo::Gemv::CompactMKL> solveblk;
          {
            Timer timer("50 SolveBlockTridiagMatrices");
            timer.reset();
            Kokkos::fence();
            for (ordinal_type i=0;i<niter;++i) {
              solveblk.run(0, T_device, x_device, b_device);
              dontopt += i;
            }          
            Kokkos::fence();
            t_solve = timer.seconds();
          }
          if (input.check) TEST_ASSERT(solveblk.check(T_org_device, b_device), success);
#elif defined(__KOKKOSBATCHED_NVIDIA_CUBLAS__)
	  std::cout << "CUBLAS compact version does not exist\n"; 
#else
	  std::cout << "TPLs (CompactMKL or CUBLAS) are not found\n"; 
#endif 
        } else {
          SolveBlockTridiagMatrices<DeviceSpace,
	    ValueType, 
	    DeviceArrayLayout,
	    VectorLength,
	    AlgoTrsv,
	    AlgoGemv> solveblk;
          {
            Timer timer("50 SolveBlockTridiagMatrices");
            timer.reset();
            Kokkos::fence();
            for (ordinal_type i=0;i<niter;++i) {
              solveblk.run(ops, T_device, x_device, b_device);
              dontopt += i;
            }          
            Kokkos::fence();
            t_solve = timer.seconds();
          }
          if (input.check) TEST_ASSERT(solveblk.check(T_org_device, b_device), success);
        }
      }
      
      const double t_matvec_per_iter = t_matvec/double(niter), t_solve_per_iter = t_solve/double(niter);
      std::cout << " matvec     = " << t_matvec_per_iter << std::endl; 
      std::cout << " extract    = " << t_extract        << " extract/matvec = " << (t_extract/t_matvec_per_iter) << std::endl; 
      //std::cout << " factor     = " << t_factorize      << " factor/matvec  = " << (t_factorize/t_matvec_per_iter) << std::endl; 
      std::cout << " factor     = " << t_factorize      << " factor/matvec  = " << (t_factorize/t_matvec_per_iter) << " flop = " << f_factorize << " flop/s = " << (f_factorize/t_factorize) << std::endl; 
      std::cout << " solve      = " << t_solve_per_iter << " solve/matvec   = " << (t_solve_per_iter/t_matvec_per_iter) << std::endl; 
      std::cout << " memory used     = " << (memsize_A + memsize_T) << std::endl; 

      return dontopt + success;
    }

  }
}
