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

#ifndef KOKKOS_BLAS3_GEMM_IMPL_HPP_
#define KOKKOS_BLAS3_GEMM_IMPL_HPP_

#include<Kokkos_Core.hpp>

namespace KokkosBlas {
namespace Impl {

// Choose Iteration Layout for copying data from global memory into scratch
// On CPUs it is more important to have consecutive write,
// On GPUs it is more important to not jump around in global memory, i.e. have coallesced loads
template<class ExecSpace, class LayoutA, class LayoutAScratch>
struct impl_gemm_choose_copy_layout {
  typedef LayoutAScratch type;
};

#ifdef KOKKOS_ENABLE_CUDA
template<class LayoutA, class LayoutAScratch>
struct impl_gemm_choose_copy_layout<Kokkos::Cuda,LayoutA,LayoutAScratch> {
  typedef LayoutA type;
};
#endif

// DeepCopy matrix block into scratch
template<class TeamHandle, class ViewTypeScratch, class ViewType, class Layout, int blockDim_i, int blockDim_j, int Transpose>
struct impl_deep_copy_matrix_block;

template<class TeamHandle, class ViewTypeScratch, class ViewType, class Layout, int blockDim_i, int blockDim_j>
struct impl_deep_copy_matrix_block<TeamHandle,ViewTypeScratch,ViewType,Layout,blockDim_i,blockDim_j,0> {
  typedef typename ViewType::non_const_value_type value_type;
  typedef Kokkos::Details::ArithTraits<value_type>     ATV;

  KOKKOS_INLINE_FUNCTION
  static void copy(const TeamHandle& team, const ViewTypeScratch& A_scr, const ViewType& A, const int& offset_i, const int& offset_j) {
    if(offset_i + blockDim_i <= A.extent_int(0) && offset_j + blockDim_j <= A.extent_int(1)) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,blockDim_j), [&] (const int j) {
        const int idx_j = offset_j+j;
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,blockDim_i), [&] (const int i) {
          const int idx_i = offset_i+i;
          A_scr(i,j) = A(idx_i,idx_j);
        });
      });
    } else {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,blockDim_j), [&] (const int j) {
        const int idx_j = offset_j+j;
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,blockDim_i), [&] (const int i) {
          const int idx_i = offset_i+i;
          A_scr(i,j) = idx_i<A.extent_int(0) && idx_j<A.extent_int(1) ? A(idx_i,idx_j) : ATV::zero();
        });
      });
    }
  }
};

template<class TeamHandle, class ViewTypeScratch, class ViewType,  int blockDim_i, int blockDim_j>
struct impl_deep_copy_matrix_block<TeamHandle,ViewTypeScratch,ViewType,Kokkos::LayoutRight,blockDim_i,blockDim_j,0> {
  typedef typename ViewType::non_const_value_type value_type;
  typedef Kokkos::Details::ArithTraits<value_type>     ATV;

  KOKKOS_INLINE_FUNCTION
  static void copy(const TeamHandle& team, const ViewTypeScratch& A_scr, const ViewType& A, const int& offset_i, const int& offset_j) {
    if(offset_i + blockDim_i <= A.extent_int(0) && offset_j + blockDim_j <= A.extent_int(1)) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,blockDim_i), [&] (const int i) {
        const int idx_i = offset_i+i;
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,blockDim_j), [&] (const int j) {
          const int idx_j = offset_j+j;
          A_scr(i,j) = A(idx_i,idx_j);
        });
      });
    } else {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,blockDim_i), [&] (const int i) {
        const int idx_i = offset_i+i;
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,blockDim_j), [&] (const int j) {
          const int idx_j = offset_j+j;
          A_scr(i,j) = idx_i<A.extent_int(0) && idx_j<A.extent_int(1) ? A(idx_i,idx_j) : ATV::zero();
        });
      });
    }
  }
};

template<class TeamHandle, class ViewTypeScratch, class ViewType, class Layout, int blockDim_i, int blockDim_j>
struct impl_deep_copy_matrix_block<TeamHandle,ViewTypeScratch,ViewType,Layout,blockDim_i,blockDim_j,1> {
  typedef typename ViewType::non_const_value_type value_type;
  typedef Kokkos::Details::ArithTraits<value_type>     ATV;

  KOKKOS_INLINE_FUNCTION
  static void copy(const TeamHandle& team, const ViewTypeScratch& A_scr, const ViewType& A, const int& offset_i, const int& offset_j) {
    if(offset_i + blockDim_i <= A.extent_int(1) && offset_j + blockDim_j <= A.extent_int(0)) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,blockDim_j), [&] (const int j) {
        const int idx_j = offset_j+j;
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,blockDim_i), [&] (const int i) {
          const int idx_i = offset_i+i;
          A_scr(i,j) = A(idx_j,idx_i);
        });
      });
    } else {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,blockDim_j), [&] (const int j) {
        const int idx_j = offset_j+j;
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,blockDim_i), [&] (const int i) {
          const int idx_i = offset_i+i;
          A_scr(i,j) = idx_i<A.extent_int(1) && idx_j<A.extent_int(0) ? A(idx_j,idx_i) : ATV::zero();
        });
      });
    }
  }
};

template<class TeamHandle, class ViewTypeScratch, class ViewType,  int blockDim_i, int blockDim_j>
struct impl_deep_copy_matrix_block<TeamHandle,ViewTypeScratch,ViewType,Kokkos::LayoutRight,blockDim_i,blockDim_j,1> {
  typedef typename ViewType::non_const_value_type value_type;
  typedef Kokkos::Details::ArithTraits<value_type>     ATV;

  KOKKOS_INLINE_FUNCTION
  static void copy(const TeamHandle& team, const ViewTypeScratch& A_scr, const ViewType& A, const int& offset_i, const int& offset_j) {
    if(offset_i + blockDim_i <= A.extent_int(1) && offset_j + blockDim_j <= A.extent_int(0)) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,blockDim_i), [&] (const int i) {
        const int idx_i = offset_i+i;
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,blockDim_j), [&] (const int j) {
          const int idx_j = offset_j+j;
          A_scr(i,j) = A(idx_j,idx_i);
        });
      });
    } else {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,blockDim_i), [&] (const int i) {
        const int idx_i = offset_i+i;
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,blockDim_j), [&] (const int j) {
          const int idx_j = offset_j+j;
          A_scr(i,j) = idx_i<A.extent_int(1) && idx_j<A.extent_int(0) ? A(idx_j,idx_i) : ATV::zero();
        });
      });
    }
  }
};

template<class TeamHandle, class ViewTypeScratch, class ViewType, class Layout, int blockDim_i, int blockDim_j>
struct impl_deep_copy_matrix_block<TeamHandle,ViewTypeScratch,ViewType,Layout,blockDim_i,blockDim_j,2> {
  typedef typename ViewType::non_const_value_type value_type;
  typedef Kokkos::Details::ArithTraits<value_type>     ATV;

  KOKKOS_INLINE_FUNCTION
  static void copy(const TeamHandle& team, const ViewTypeScratch& A_scr, const ViewType& A, const int& offset_i, const int& offset_j) {
    if(offset_i + blockDim_i <= A.extent_int(1) && offset_j + blockDim_j <= A.extent_int(0)) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,blockDim_j), [&] (const int j) {
        const int idx_j = offset_j+j;
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,blockDim_i), [&] (const int i) {
          const int idx_i = offset_i+i;
          A_scr(i,j) = ATV::conj(A(idx_j,idx_i));
        });
      });
    } else {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,blockDim_j), [&] (const int j) {
        const int idx_j = offset_j+j;
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,blockDim_i), [&] (const int i) {
          const int idx_i = offset_i+i;
          A_scr(i,j) = idx_i<A.extent_int(1) && idx_j<A.extent_int(0) ? ATV::conj(A(idx_j,idx_i)) : ATV::zero();
        });
      });
    }
  }
};

template<class TeamHandle, class ViewTypeScratch, class ViewType,  int blockDim_i, int blockDim_j>
struct impl_deep_copy_matrix_block<TeamHandle,ViewTypeScratch,ViewType,Kokkos::LayoutRight,blockDim_i,blockDim_j,2> {
  typedef typename ViewType::non_const_value_type value_type;
  typedef Kokkos::Details::ArithTraits<value_type>     ATV;

  KOKKOS_INLINE_FUNCTION
  static void copy(const TeamHandle& team, const ViewTypeScratch& A_scr, const ViewType& A, const int& offset_i, const int& offset_j) {
    if(offset_i + blockDim_i <= A.extent_int(1) && offset_j + blockDim_j <= A.extent_int(0)) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,blockDim_i), [&] (const int i) {
        const int idx_i = offset_i+i;
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,blockDim_j), [&] (const int j) {
          const int idx_j = offset_j+j;
          A_scr(i,j) = A(idx_j,idx_i);
        });
      });
    } else {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,blockDim_i), [&] (const int i) {
        const int idx_i = offset_i+i;
        Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,blockDim_j), [&] (const int j) {
          const int idx_j = offset_j+j;
          A_scr(i,j) = idx_i<A.extent_int(1) && idx_j<A.extent_int(0) ? ATV::conj(A(idx_j,idx_i)) : ATV::zero();
        });
      });
    }
  }
};

template<class TeamHandle, class ViewType, class ViewTypeScratch, class Layout, int blockDim_i, int blockDim_j>
struct impl_update_matrix_block {
  typedef typename ViewType::non_const_value_type value_type;
  typedef Kokkos::Details::ArithTraits<value_type>     ATV;

  KOKKOS_INLINE_FUNCTION
  static void update(const TeamHandle& team, const value_type& beta , const ViewType& A,
                                             const value_type& alpha, const ViewTypeScratch& A_scr,
                                             const int& offset_i, const int& offset_j) {
    if(offset_i + blockDim_i <= A.extent_int(0) && offset_j + blockDim_j <= A.extent_int(1)) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,blockDim_j), [&] (const int j) {
        const int idx_j = offset_j+j;
        if(beta == ATV::zero()) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,blockDim_i), [&] (const int i) {
            const int idx_i = offset_i+i;
            A(idx_i,idx_j) = alpha * A_scr(i,j);
          });
        } else {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,blockDim_i), [&] (const int i) {
            const int idx_i = offset_i+i;
            A(idx_i,idx_j) = beta * A(idx_i,idx_j) + alpha * A_scr(i,j);
          });
        }
      });
    } else {
      const int range_i = offset_i + blockDim_i <= A.extent_int(0)?blockDim_i:A.extent_int(0)%blockDim_i;
      const int range_j = offset_j + blockDim_j <= A.extent_int(1)?blockDim_j:A.extent_int(1)%blockDim_j;
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,range_j), [&] (const int j) {
        const int idx_j = offset_j+j;
        if(beta == ATV::zero()) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,range_i), [&] (const int i) {
            const int idx_i = offset_i+i;
            A(idx_i,idx_j) = alpha * A_scr(i,j);
          });
        } else {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,range_i), [&] (const int i) {
            const int idx_i = offset_i+i;
            A(idx_i,idx_j) = beta * A(idx_i,idx_j) + alpha * A_scr(i,j);
          });
        }
      });
    }
  }
};

template<class TeamHandle, class ViewType, class ViewTypeScratch, int blockDim_i, int blockDim_j>
struct impl_update_matrix_block<TeamHandle,ViewType,ViewTypeScratch,Kokkos::LayoutRight,blockDim_i,blockDim_j> {
  typedef typename ViewType::non_const_value_type value_type;
  typedef Kokkos::Details::ArithTraits<value_type>     ATV;

  KOKKOS_INLINE_FUNCTION
  static void update(const TeamHandle& team, const value_type& beta , const ViewType& A,
                                             const value_type& alpha, const ViewTypeScratch& A_scr,
                                             const int& offset_i, const int& offset_j) {
    if(offset_i + blockDim_i <= A.extent_int(0) && offset_j + blockDim_j <= A.extent_int(1)) {
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,blockDim_i), [&] (const int i) {
        const int idx_i = offset_i+i;
        if(beta == ATV::zero()) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,blockDim_j), [&] (const int j) {
            const int idx_j = offset_j+j;
            A(idx_i,idx_j) = alpha * A_scr(i,j);
          });
        } else {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,blockDim_j), [&] (const int j) {
            const int idx_j = offset_j+j;
            A(idx_i,idx_j) = beta * A(idx_i,idx_j) + alpha * A_scr(i,j);
          });
        }
      });
    } else {
      const int range_i = offset_i + blockDim_i <= A.extent_int(0)?blockDim_i:A.extent_int(0)%blockDim_i;
      const int range_j = offset_j + blockDim_j <= A.extent_int(1)?blockDim_j:A.extent_int(1)%blockDim_j;
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team,range_i), [&] (const int i) {
        const int idx_i = offset_i+i;
        if(beta == ATV::zero()) {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,range_j), [&] (const int j) {
            const int idx_j = offset_j+j;
            A(idx_i,idx_j) = alpha * A_scr(i,j);
          });
        } else {
          Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,range_j), [&] (const int j) {
            const int idx_j = offset_j+j;
            A(idx_i,idx_j) = beta * A(idx_i,idx_j) + alpha * A_scr(i,j);
          });

        }
      });
    }
  }
};

// Compute a single A block 8 B block, also do an in-place no-additional blocking team GEMM
template<class TeamHandle, class ViewTypeA, class ViewTypeB, class ViewTypeC>
KOKKOS_INLINE_FUNCTION
void impl_team_gemm_block(const TeamHandle& team, const ViewTypeC& C, const ViewTypeA& A, const ViewTypeB& B) {
  typedef typename ViewTypeC::non_const_value_type ScalarC;
// GNU COMPILER BUG WORKAROUND
#if defined(KOKKOS_COMPILER_GNU) || !defined(__CUDA_ARCH__)
  int blockA0 = A.extent_int(0);
  int blockA1 = A.extent_int(1);
  int blockB1 = B.extent_int(1);
#else
  const int blockA0 = A.extent_int(0);
  const int blockA1 = A.extent_int(1);
  const int blockB1 = B.extent_int(1);
#endif
  Kokkos::parallel_for(Kokkos::TeamThreadRange(team,blockA0), [&] (const int i) {
#if defined(__CUDA_ARCH__) || !defined(KOKKOS_ENABLE_OPENMP)
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,blockB1/4), [&] (const int B_j) {
#else
  #if defined(KOKKOS_COMPILER_GNU)
    #if (KOKKOS_COMPILER_GNU > 485 )
    #pragma omp simd
    #endif
  #else
    #pragma omp simd
  #endif
    for(int B_j=0; B_j<blockB1/4; B_j++) {
#endif
      ScalarC C_ij0 = 0;
      ScalarC C_ij1 = 0;
      ScalarC C_ij2 = 0;
      ScalarC C_ij3 = 0;
      for(int j = 0; j < blockA1; j++) {
        ScalarC A_ij = A(i,j);
        C_ij0 += A_ij*B(j,B_j);
        C_ij1 += A_ij*B(j,B_j+blockB1/4);
        C_ij2 += A_ij*B(j,B_j+2*blockB1/4);
        C_ij3 += A_ij*B(j,B_j+3*blockB1/4);
      }
      C(i,B_j) += C_ij0;
      C(i,B_j+blockB1/4) += C_ij1;
      C(i,B_j+2*blockB1/4) += C_ij2;
      C(i,B_j+3*blockB1/4) += C_ij3;
#if defined(__CUDA_ARCH__) || !defined(KOKKOS_ENABLE_OPENMP)
    });
#else
    }
#endif
  });
}

template<int TransposeA, int TransposeB>
struct impl_gemm_label;

template<>
struct impl_gemm_label<0,0> {
  static constexpr const char* label = "KokkosBlas::gemm[NN]";
};
template<>
struct impl_gemm_label<0,1> {
  static constexpr const char* label = "KokkosBlas::gemm[NT]";
};
template<>
struct impl_gemm_label<0,2> {
  static constexpr const char* label = "KokkosBlas::gemm[NC]";
};

template<>
struct impl_gemm_label<1,0> {
  static constexpr const char* label = "KokkosBlas::gemm[TN]";
};
template<>
struct impl_gemm_label<1,1> {
  static constexpr const char* label = "KokkosBlas::gemm[TT]";
};
template<>
struct impl_gemm_label<1,2> {
  static constexpr const char* label = "KokkosBlas::gemm[TC]";
};

template<>
struct impl_gemm_label<2,0> {
  static constexpr const char* label = "KokkosBlas::gemm[CN]";
};
template<>
struct impl_gemm_label<2,1> {
  static constexpr const char* label = "KokkosBlas::gemm[CT]";
};
template<>
struct impl_gemm_label<2,2> {
  static constexpr const char* label = "KokkosBlas::gemm[CC]";
};


template<class ExecSpace, class ViewTypeA, class ViewTypeB, class ViewTypeC,
          int blockA0, int blockA1, int blockB1, int TransposeA, int TransposeB>
struct GEMMImpl {
  ViewTypeA A;
  ViewTypeB B;
  ViewTypeC C;
  typedef typename ViewTypeA::non_const_value_type ScalarA;
  typedef typename ViewTypeB::non_const_value_type ScalarB;
  typedef typename ViewTypeC::non_const_value_type ScalarC;

  const int num_blocks_0;
  const int num_blocks_1;
  int scratch_level;

  ScalarC alpha, beta;
  typedef Kokkos::View<ScalarA[blockA0][blockA1],Kokkos::LayoutLeft,typename ExecSpace::scratch_memory_space>
    ViewTypeAScratch;
  typedef Kokkos::View<ScalarB[blockA1][blockB1],Kokkos::LayoutRight,typename ExecSpace::scratch_memory_space>
    ViewTypeBScratch;
  typedef Kokkos::View<ScalarC[blockA0][blockB1],Kokkos::LayoutRight,typename ExecSpace::scratch_memory_space>
    ViewTypeCScratch;

  GEMMImpl(const ScalarA& alpha_, const ViewTypeA& A_, const ViewTypeB& B_, const ScalarC& beta_, const ViewTypeC& C_):A(A_),B(B_),C(C_),
      num_blocks_0((C.extent_int(0)+blockA0-1)/blockA0),num_blocks_1((C.extent_int(1)+blockB1-1)/blockB1) {
    scratch_level = 0;
    alpha = alpha_;
    beta = beta_;
  }

  void run(int team_size, int vector_length, int scr_level) {
    scratch_level = scr_level;
    int scratch_memory_size =
      ViewTypeAScratch::shmem_size() +
      ViewTypeBScratch::shmem_size() +
      ViewTypeCScratch::shmem_size();

    Kokkos::TeamPolicy<ExecSpace,Kokkos::LaunchBounds<384,2>> policy(num_blocks_0*num_blocks_1,team_size,vector_length);

    Kokkos::parallel_for(impl_gemm_label<TransposeA,TransposeB>::label,policy.set_scratch_size(scratch_level,Kokkos::PerTeam(scratch_memory_size)),*this);
  }

  KOKKOS_INLINE_FUNCTION
  void operator() (const typename Kokkos::TeamPolicy<ExecSpace>::member_type& team) const {
    // This team is responsible for computing a single block of C
    const int league_rank = team.league_rank();
    const int num_blocks = num_blocks_1;
    const int i_offset = (league_rank/num_blocks)*blockA0;
    const int j_offset = (league_rank%num_blocks)*blockB1;

    ViewTypeAScratch A_scr(team.team_scratch(scratch_level));
    ViewTypeBScratch B_scr(team.team_scratch(scratch_level));
    ViewTypeCScratch C_scr(team.team_scratch(scratch_level));
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team,blockA0), [&] (const int i) {
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(team,blockB1), [&] (const int j) {
        C_scr(i,j) = 0;
      });
    });
    team.team_barrier();

    // Move along the inner dimension in blocks
    const int length = TransposeA>0?A.extent_int(0):A.extent_int(1);
    for(int A_j = 0; A_j < length; A_j += blockA1) {
      // Load A block into scratch

      impl_deep_copy_matrix_block<typename Kokkos::TeamPolicy<ExecSpace>::member_type,
                                  ViewTypeAScratch,ViewTypeA,
                                  typename impl_gemm_choose_copy_layout<ExecSpace,
                                     typename ViewTypeA::array_layout,
                                     typename ViewTypeAScratch::array_layout>::type,
                                  blockA0,blockA1,TransposeA>::copy(team,A_scr,A,i_offset,A_j);

      // Load B block into scratch
      impl_deep_copy_matrix_block<typename Kokkos::TeamPolicy<ExecSpace>::member_type,
                                  ViewTypeBScratch,ViewTypeB,
                                  typename impl_gemm_choose_copy_layout<ExecSpace,
                                     typename ViewTypeB::array_layout,
                                     typename ViewTypeBScratch::array_layout>::type,
                                  blockA1,blockB1,TransposeB>::copy(team,B_scr,B,A_j,j_offset);

      // Wait for A and B block to be in scratch memory
      team.team_barrier();

      // Add contribution from multiplying the A and B block to the C block
      impl_team_gemm_block(team,C_scr,A_scr,B_scr);

      // Wait for subblock computation to be done before loading the next A and B block
      team.team_barrier();
    }
    // Write back the C block from scratch to main memory
    impl_update_matrix_block<typename Kokkos::TeamPolicy<ExecSpace>::member_type,
                                      ViewTypeC,ViewTypeCScratch,
                                      typename ViewTypeC::array_layout,
                                      blockA0,blockB1>::update(team,beta,C,alpha,C_scr,i_offset,j_offset);
  }
};

}
}
#endif
