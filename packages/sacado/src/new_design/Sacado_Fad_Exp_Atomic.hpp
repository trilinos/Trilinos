// @HEADER
// ***********************************************************************
//
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
// USA
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
//
// ***********************************************************************
// @HEADER

#ifndef SACADO_FAD_EXP_ATOMIC_HPP
#define SACADO_FAD_EXP_ATOMIC_HPP

#include "Sacado_ConfigDefs.h"
#if defined(HAVE_SACADO_KOKKOSCORE)

#include "Sacado_Fad_Exp_ViewFad.hpp"
#include "Kokkos_Atomic.hpp"
#include "impl/Kokkos_Error.hpp"

namespace Sacado {

  namespace Fad {
  namespace Exp {

    // Overload of Kokkos::atomic_add for ViewFad types.
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION
    void atomic_add(ViewFadPtr<ValT,sl,ss,U> dst, const Expr<T>& xx) {
      using Kokkos::atomic_add;

      const typename Expr<T>::derived_type& x = xx.derived();

      const int xsz = x.size();
      const int sz = dst->size();

      // We currently cannot handle resizing since that would need to be
      // done atomically.
      if (xsz > sz)
        Kokkos::abort(
          "Sacado error: Fad resize within atomic_add() not supported!");

      if (xsz != sz && sz > 0 && xsz > 0)
        Kokkos::abort(
          "Sacado error: Fad assignment of incompatiable sizes!");


      if (sz > 0 && xsz > 0) {
        SACADO_FAD_DERIV_LOOP(i,sz)
          atomic_add(&(dst->fastAccessDx(i)), x.fastAccessDx(i));
      }
      SACADO_FAD_THREAD_SINGLE
        atomic_add(&(dst->val()), x.val());
    }

    namespace Impl {
      // Our implementation of Kokkos::atomic_oper_fetch() and
      // Kokkos::atomic_fetch_oper() for Sacado types
      template <typename Oper, typename DestPtrT, typename ValT, typename T>
      SACADO_INLINE_FUNCTION
      typename Sacado::BaseExprType< Expr<T> >::type
      atomic_oper_fetch_impl(const Oper& op, DestPtrT dest, ValT* dest_val,
                             const Expr<T>& x)
      {
        typedef typename Sacado::BaseExprType< Expr<T> >::type return_type;
        const typename Expr<T>::derived_type& val = x.derived();

#if defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST)
        while (!Kokkos::Impl::lock_address_host_space((void*)dest_val))
          ;
        Kokkos::memory_fence();
        return_type return_val = op.apply(*dest, val);
        *dest                  = return_val;
        Kokkos::memory_fence();
        Kokkos::Impl::unlock_address_host_space((void*)dest_val);
        return return_val;
#elif defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA)
	// It is not allowed to define SACADO_VIEW_CUDA_HIERARCHICAL or
	// SACADO_VIEW_CUDA_HIERARCHICAL_DFAD and use Sacado inside a team-based
	// kernel without Sacado hierarchical parallelism.  So use the
	// team-based version only if blockDim.x > 1 (i.e., a team policy)
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) || defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
	const bool use_team = (blockDim.x > 1);
#else
	const bool use_team = false;
#endif
	if (use_team) {
	  int go = 1;
	  while (go) {
	    if (threadIdx.x == 0)
	      go = !Kokkos::Impl::lock_address_cuda_space((void*)dest_val);
	    go = Kokkos::shfl(go, 0, blockDim.x);
	  }
	  Kokkos::memory_fence();
	  return_type return_val = op.apply(*dest, val);
	  *dest                  = return_val;
	  Kokkos::memory_fence();
	  if (threadIdx.x == 0)
	    Kokkos::Impl::unlock_address_cuda_space((void*)dest_val);
	  return return_val;
	}
	else {
	  return_type return_val;
	  // This is a way to (hopefully) avoid dead lock in a warp
	  int done                 = 0;
	  unsigned int mask        =  __activemask() ;
	  unsigned int active      = __ballot_sync(mask, 1);
	  unsigned int done_active = 0;
	  while (active != done_active) {
	    if (!done) {
	      if (Kokkos::Impl::lock_address_cuda_space((void*)dest_val)) {
		Kokkos::memory_fence();
		return_val = op.apply(*dest, val);
		*dest      = return_val;
		Kokkos::memory_fence();
		Kokkos::Impl::unlock_address_cuda_space((void*)dest_val);
		done = 1;
	      }
	    }
	    done_active = __ballot_sync(mask, done);
	  }
	  return return_val;
	}
#elif defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HIP_GPU)
	// It is not allowed to define SACADO_VIEW_CUDA_HIERARCHICAL or
	// SACADO_VIEW_CUDA_HIERARCHICAL_DFAD and use Sacado inside a team-based
	// kernel without Sacado hierarchical parallelism.  So use the
	// team-based version only if blockDim.x > 1 (i.e., a team policy)
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) || defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
	const bool use_team = (blockDim.x > 1);
#else
	const bool use_team = false;
#endif
	if (use_team) {
	  int go = 1;
	  while (go) {
	    if (threadIdx.x == 0)
	      go = !Kokkos::Impl::lock_address_hip_space((void*)dest_val);
	    go = Kokkos::Experimental::shfl(go, 0, blockDim.x);
	  }
	  Kokkos::memory_fence();
	  return_type return_val = op.apply(*dest, val);
	  *dest                  = return_val;
	  Kokkos::memory_fence();
	  if (threadIdx.x == 0)
	    Kokkos::Impl::unlock_address_hip_space((void*)dest_val);
	  return return_val;
	}
	else {
	  return_type return_val;
	  int done                 = 0;
	  unsigned int active      = __ballot(1);
	  unsigned int done_active = 0;
	  while (active != done_active) {
	    if (!done) {
	      if (Kokkos::Impl::lock_address_hip_space((void*)dest_val)) {
		return_val = op.apply(*dest, val);
		*dest      = return_val;
		Kokkos::Impl::unlock_address_hip_space((void*)dest_val);
		done = 1;
	      }
	    }
	    done_active = __ballot(done);
	  }
	  return return_val;
	}
#endif
      }

      template <typename Oper, typename DestPtrT, typename ValT, typename T>
      SACADO_INLINE_FUNCTION
      typename Sacado::BaseExprType< Expr<T> >::type
      atomic_fetch_oper_impl(const Oper& op, DestPtrT dest, ValT* dest_val,
                             const Expr<T>& x)
      {
        typedef typename Sacado::BaseExprType< Expr<T> >::type return_type;
        const typename Expr<T>::derived_type& val = x.derived();

#ifdef KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HOST
        while (!Kokkos::Impl::lock_address_host_space((void*)dest_val))
          ;
        Kokkos::memory_fence();
        return_type return_val = *dest;
        *dest                  = op.apply(return_val, val);
        Kokkos::memory_fence();
        Kokkos::Impl::unlock_address_host_space((void*)dest_val);
        return return_val;
#elif defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_CUDA)
	// It is not allowed to define SACADO_VIEW_CUDA_HIERARCHICAL or
	// SACADO_VIEW_CUDA_HIERARCHICAL_DFAD and use Sacado inside a team-based
	// kernel without Sacado hierarchical parallelism.  So use the
	// team-based version only if blockDim.x > 1 (i.e., a team policy)
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) || defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
	const bool use_team = (blockDim.x > 1);
#else
	const bool use_team = false;
#endif
	if (use_team) {
	  int go = 1;
	  while (go) {
	    if (threadIdx.x == 0)
	      go = !Kokkos::Impl::lock_address_cuda_space((void*)dest_val);
	    go = Kokkos::shfl(go, 0, blockDim.x);
	  }
	  Kokkos::memory_fence();
	  return_type return_val = *dest;
	  *dest                  = op.apply(return_val, val);
	  Kokkos::memory_fence();
	  if (threadIdx.x == 0)
	    Kokkos::Impl::unlock_address_cuda_space((void*)dest_val);
	  return return_val;
	}
	else {
	  return_type return_val;
	  // This is a way to (hopefully) avoid dead lock in a warp
	  int done                 = 0;
	  unsigned int mask        =  __activemask() ;
	  unsigned int active      = __ballot_sync(mask, 1);
	  unsigned int done_active = 0;
	  while (active != done_active) {
	    if (!done) {
	      if (Kokkos::Impl::lock_address_cuda_space((void*)dest_val)) {
		Kokkos::memory_fence();
		return_val = *dest;
		*dest      = op.apply(return_val, val);
		Kokkos::memory_fence();
		Kokkos::Impl::unlock_address_cuda_space((void*)dest_val);
		done = 1;
	      }
	    }
	    done_active = __ballot_sync(mask, done);
	  }
	  return return_val;
	}
#elif defined(KOKKOS_ACTIVE_EXECUTION_MEMORY_SPACE_HIP_GPU)
	// It is not allowed to define SACADO_VIEW_CUDA_HIERARCHICAL or
	// SACADO_VIEW_CUDA_HIERARCHICAL_DFAD and use Sacado inside a team-based
	// kernel without Sacado hierarchical parallelism.  So use the
	// team-based version only if blockDim.x > 1 (i.e., a team policy)
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) || defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
	const bool use_team = (blockDim.x > 1);
#else
	const bool use_team = false;
#endif
	if (use_team) {
	  int go = 1;
	  while (go) {
	    if (threadIdx.x == 0)
	      go = !Kokkos::Impl::lock_address_hip_space((void*)dest_val);
	    go = Kokkos::Experimental::shfl(go, 0, blockDim.x);
	  }
	  Kokkos::memory_fence();
	  return_type return_val = *dest;
	  *dest                  = op.apply(return_val, val);
	  Kokkos::memory_fence();
	  if (threadIdx.x == 0)
	    Kokkos::Impl::unlock_address_hip_space((void*)dest_val);
	  return return_val;
	}
	else {
	  return_type return_val;
	  int done                 = 0;
	  unsigned int active      = __ballot(1);
	  unsigned int done_active = 0;
	  while (active != done_active) {
	    if (!done) {
	      if (Kokkos::Impl::lock_address_hip_space((void*)dest_val)) {
		return_val = *dest;
		*dest      = op.apply(return_val, val);
		Kokkos::Impl::unlock_address_hip_space((void*)dest_val);
		done = 1;
	      }
	    }
	    done_active = __ballot(done);
	  }
	  return return_val;
	}
#endif
      }

      // Overloads of Kokkos::atomic_oper_fetch/Kokkos::atomic_fetch_oper
      // for Sacado types
      template <typename Oper, typename S>
      SACADO_INLINE_FUNCTION GeneralFad<S>
      atomic_oper_fetch(const Oper& op, GeneralFad<S>* dest,
                        const GeneralFad<S>& val)
      {
        return Impl::atomic_oper_fetch_impl(op, dest, &(dest->val()), val);
      }
      template <typename Oper, typename ValT, unsigned sl, unsigned ss,
                typename U, typename T>
      SACADO_INLINE_FUNCTION U
      atomic_oper_fetch(const Oper& op, ViewFadPtr<ValT,sl,ss,U> dest,
                        const Expr<T>& val)
      {
        return Impl::atomic_oper_fetch_impl(op, dest, &dest.val(), val);
      }

      template <typename Oper, typename S>
      SACADO_INLINE_FUNCTION GeneralFad<S>
      atomic_fetch_oper(const Oper& op, GeneralFad<S>* dest,
                        const GeneralFad<S>& val)
      {
        return Impl::atomic_fetch_oper_impl(op, dest, &(dest->val()), val);
      }
      template <typename Oper, typename ValT, unsigned sl, unsigned ss,
                typename U, typename T>
      SACADO_INLINE_FUNCTION U
      atomic_fetch_oper(const Oper& op, ViewFadPtr<ValT,sl,ss,U> dest,
                        const Expr<T>& val)
      {
        return Impl::atomic_fetch_oper_impl(op, dest, &dest.val(), val);
      }

      // Our definition of the various Oper classes to be more type-flexible
      struct MaxOper {
        template <class Scalar1, class Scalar2>
        KOKKOS_FORCEINLINE_FUNCTION
        static auto apply(const Scalar1& val1, const Scalar2& val2)
          -> decltype(max(val1,val2))
        {
          return max(val1,val2);
        }
      };
      struct MinOper {
        template <class Scalar1, class Scalar2>
        KOKKOS_FORCEINLINE_FUNCTION
        static auto apply(const Scalar1& val1, const Scalar2& val2)
          -> decltype(min(val1,val2))
        {
          return min(val1,val2);
        }
      };
      struct AddOper {
        template <class Scalar1, class Scalar2>
        KOKKOS_FORCEINLINE_FUNCTION
        static auto apply(const Scalar1& val1, const Scalar2& val2)
          -> decltype(val1+val2)
        {
          return val1 + val2;
        }
      };
      struct SubOper {
        template <class Scalar1, class Scalar2>
        KOKKOS_FORCEINLINE_FUNCTION
        static auto apply(const Scalar1& val1, const Scalar2& val2)
          -> decltype(val1-val2)
        {
          return val1 - val2;
        }
      };
      struct MulOper {
        template <class Scalar1, class Scalar2>
        KOKKOS_FORCEINLINE_FUNCTION
        static auto apply(const Scalar1& val1, const Scalar2& val2)
          -> decltype(val1*val2)
        {
          return val1 * val2;
        }
      };
      struct DivOper {
        template <class Scalar1, class Scalar2>
        KOKKOS_FORCEINLINE_FUNCTION
        static auto apply(const Scalar1& val1, const Scalar2& val2)
          -> decltype(val1/val2)
        {
          return val1 / val2;
        }
      };

    } // Impl

    // Overload of Kokkos::atomic_*_fetch() and Kokkos::atomic_fetch_*()
    // for Sacado types
    template <typename S>
    SACADO_INLINE_FUNCTION GeneralFad<S>
    atomic_max_fetch(GeneralFad<S>* dest, const GeneralFad<S>& val) {
      return Impl::atomic_oper_fetch(Impl::MaxOper(), dest, val);
    }
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION U
    atomic_max_fetch(ViewFadPtr<ValT,sl,ss,U> dest, const Expr<T>& val) {
      return Impl::atomic_oper_fetch(Impl::MaxOper(), dest, val);
    }
    template <typename S>
    SACADO_INLINE_FUNCTION GeneralFad<S>
    atomic_min_fetch(GeneralFad<S>* dest, const GeneralFad<S>& val) {
      return Impl::atomic_oper_fetch(Impl::MinOper(), dest, val);
    }
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION U
    atomic_min_fetch(ViewFadPtr<ValT,sl,ss,U> dest, const Expr<T>& val) {
      return Impl::atomic_oper_fetch(Impl::MinOper(), dest, val);
    }
    template <typename S>
    SACADO_INLINE_FUNCTION GeneralFad<S>
    atomic_add_fetch(GeneralFad<S>* dest, const GeneralFad<S>& val) {
      return Impl::atomic_oper_fetch(Impl::AddOper(), dest, val);
    }
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION U
    atomic_add_fetch(ViewFadPtr<ValT,sl,ss,U> dest, const Expr<T>& val) {
      return Impl::atomic_oper_fetch(Impl::AddOper(), dest, val);
    }
    template <typename S>
    SACADO_INLINE_FUNCTION GeneralFad<S>
    atomic_sub_fetch(GeneralFad<S>* dest, const GeneralFad<S>& val) {
      return Impl::atomic_oper_fetch(Impl::SubOper(), dest, val);
    }
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION U
    atomic_sub_fetch(ViewFadPtr<ValT,sl,ss,U> dest, const Expr<T>& val) {
      return Impl::atomic_oper_fetch(Impl::SubOper(), dest, val);
    }
    template <typename S>
    SACADO_INLINE_FUNCTION GeneralFad<S>
    atomic_mul_fetch(GeneralFad<S>* dest, const GeneralFad<S>& val) {
      return atomic_oper_fetch(Impl::MulOper(), dest, val);
    }
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION U
    atomic_mul_fetch(ViewFadPtr<ValT,sl,ss,U> dest, const Expr<T>& val) {
      return Impl::atomic_oper_fetch(Impl::MulOper(), dest, val);
    }
    template <typename S>
    SACADO_INLINE_FUNCTION GeneralFad<S>
    atomic_div_fetch(GeneralFad<S>* dest, const GeneralFad<S>& val) {
      return Impl::atomic_oper_fetch(Impl::DivOper(), dest, val);
    }
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION U
    atomic_div_fetch(ViewFadPtr<ValT,sl,ss,U> dest, const Expr<T>& val) {
      return Impl::atomic_oper_fetch(Impl::DivOper(), dest, val);
    }

    template <typename S>
    SACADO_INLINE_FUNCTION GeneralFad<S>
    atomic_fetch_max(GeneralFad<S>* dest, const GeneralFad<S>& val) {
      return Impl::atomic_fetch_oper(Impl::MaxOper(), dest, val);
    }
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION U
    atomic_fetch_max(ViewFadPtr<ValT,sl,ss,U> dest, const Expr<T>& val) {
      return Impl::atomic_fetch_oper(Impl::MaxOper(), dest, val);
    }
    template <typename S>
    SACADO_INLINE_FUNCTION GeneralFad<S>
    atomic_fetch_min(GeneralFad<S>* dest, const GeneralFad<S>& val) {
      return Impl::atomic_fetch_oper(Impl::MinOper(), dest, val);
    }
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION U
    atomic_fetch_min(ViewFadPtr<ValT,sl,ss,U> dest, const Expr<T>& val) {
      return Impl::atomic_fetch_oper(Impl::MinOper(), dest, val);
    }
    template <typename S>
    SACADO_INLINE_FUNCTION GeneralFad<S>
    atomic_fetch_add(GeneralFad<S>* dest, const GeneralFad<S>& val) {
      return Impl::atomic_fetch_oper(Impl::AddOper(), dest, val);
    }
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION U
    atomic_fetch_add(ViewFadPtr<ValT,sl,ss,U> dest, const Expr<T>& val) {
      return Impl::atomic_fetch_oper(Impl::AddOper(), dest, val);
    }
    template <typename S>
    SACADO_INLINE_FUNCTION GeneralFad<S>
    atomic_fetch_sub(GeneralFad<S>* dest, const GeneralFad<S>& val) {
      return Impl::atomic_fetch_oper(Impl::SubOper(), dest, val);
    }
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION U
    atomic_fetch_sub(ViewFadPtr<ValT,sl,ss,U> dest, const Expr<T>& val) {
      return Impl::atomic_fetch_oper(Impl::SubOper(), dest, val);
    }
    template <typename S>
    SACADO_INLINE_FUNCTION GeneralFad<S>
    atomic_fetch_mul(GeneralFad<S>* dest, const GeneralFad<S>& val) {
      return Impl::atomic_fetch_oper(Impl::MulOper(), dest, val);
    }
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION U
    atomic_fetch_mul(ViewFadPtr<ValT,sl,ss,U> dest, const Expr<T>& val) {
      return Impl::atomic_fetch_oper(Impl::MulOper(), dest, val);
    }
    template <typename S>
    SACADO_INLINE_FUNCTION GeneralFad<S>
    atomic_fetch_div(GeneralFad<S>* dest, const GeneralFad<S>& val) {
      return Impl::atomic_fetch_oper(Impl::DivOper(), dest, val);
    }
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    SACADO_INLINE_FUNCTION U
    atomic_fetch_div(ViewFadPtr<ValT,sl,ss,U> dest, const Expr<T>& val) {
      return Impl::atomic_fetch_oper(Impl::DivOper(), dest, val);
    }

  } // namespace Exp
  } // namespace Fad

} // namespace Sacado

#endif // HAVE_SACADO_KOKKOSCORE
#endif // SACADO_FAD_EXP_VIEWFAD_HPP
