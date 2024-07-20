// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_EXP_ATOMIC_HPP
#define SACADO_FAD_EXP_ATOMIC_HPP

#include "Sacado_ConfigDefs.h"
#if defined(HAVE_SACADO_KOKKOS)

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
      // Kokkos::atomic_fetch_oper() for Sacado types on host
      template <typename Oper, typename DestPtrT, typename ValT, typename T>
      typename Sacado::BaseExprType< Expr<T> >::type
      atomic_oper_fetch_host(const Oper& op, DestPtrT dest, ValT* dest_val,
                             const Expr<T>& x)
      {
        typedef typename Sacado::BaseExprType< Expr<T> >::type return_type;
        const typename Expr<T>::derived_type& val = x.derived();

#ifdef KOKKOS_INTERNAL_NOT_PARALLEL
        auto scope = desul::MemoryScopeCaller();
#else
        auto scope = desul::MemoryScopeDevice();
#endif

        while (!desul::Impl::lock_address((void*)dest_val, scope))
          ;
        desul::atomic_thread_fence(desul::MemoryOrderAcquire(), scope);
        return_type return_val = op.apply(*dest, val);
        *dest                  = return_val;
        desul::atomic_thread_fence(desul::MemoryOrderRelease(), scope);
        desul::Impl::unlock_address((void*)dest_val, scope);
        return return_val;
      }

      template <typename Oper, typename DestPtrT, typename ValT, typename T>
      typename Sacado::BaseExprType< Expr<T> >::type
      atomic_fetch_oper_host(const Oper& op, DestPtrT dest, ValT* dest_val,
                             const Expr<T>& x)
      {
        typedef typename Sacado::BaseExprType< Expr<T> >::type return_type;
        const typename Expr<T>::derived_type& val = x.derived();

#ifdef KOKKOS_INTERNAL_NOT_PARALLEL
        auto scope = desul::MemoryScopeCaller();
#else
        auto scope = desul::MemoryScopeDevice();
#endif

        while (!desul::Impl::lock_address((void*)dest_val, scope))
          ;
        desul::atomic_thread_fence(desul::MemoryOrderAcquire(), scope);
        return_type return_val = *dest;
        *dest                  = op.apply(return_val, val);
        desul::atomic_thread_fence(desul::MemoryOrderRelease(), scope);
        desul::Impl::unlock_address((void*)dest_val, scope);
        return return_val;
      }

      // Helper function to decide if we are using team-based parallelism
#if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP)
      __device__
      inline bool atomics_use_team() {
#if defined(SACADO_VIEW_CUDA_HIERARCHICAL) || defined(SACADO_VIEW_CUDA_HIERARCHICAL_DFAD)
        // It is not allowed to define SACADO_VIEW_CUDA_HIERARCHICAL or
        // SACADO_VIEW_CUDA_HIERARCHICAL_DFAD and use Sacado inside a team-based
        // kernel without Sacado hierarchical parallelism.  So use the
        // team-based version only if blockDim.x > 1 (i.e., a team policy)
        return (blockDim.x > 1);
#else
        return false;
#endif
      }
#endif

#if defined(KOKKOS_ENABLE_CUDA)

      // Our implementation of Kokkos::atomic_oper_fetch() and
      // Kokkos::atomic_fetch_oper() for Sacado types on device
      template <typename Oper, typename DestPtrT, typename ValT, typename T>
      __device__
      typename Sacado::BaseExprType< Expr<T> >::type
      atomic_oper_fetch_device(const Oper& op, DestPtrT dest, ValT* dest_val,
                               const Expr<T>& x)
      {
        typedef typename Sacado::BaseExprType< Expr<T> >::type return_type;
        const typename Expr<T>::derived_type& val = x.derived();

        auto scope = desul::MemoryScopeDevice();

        if (atomics_use_team()) {
          int go = 1;
          while (go) {
            if (threadIdx.x == 0)
              go = !desul::Impl::lock_address_cuda((void*)dest_val, scope);
            go = Kokkos::shfl(go, 0, blockDim.x);
          }
          desul::atomic_thread_fence(desul::MemoryOrderAcquire(), scope);
          return_type return_val = op.apply(*dest, val);
          *dest                  = return_val;
          desul::atomic_thread_fence(desul::MemoryOrderRelease(), scope);
          if (threadIdx.x == 0)
            desul::Impl::unlock_address_cuda((void*)dest_val, scope);
          return return_val;
        }
        else {
          return_type return_val;
          // This is a way to avoid dead lock in a warp
          int done                 = 0;
          unsigned int mask        =  __activemask() ;
          unsigned int active      = __ballot_sync(mask, 1);
          unsigned int done_active = 0;
          while (active != done_active) {
            if (!done) {
              if (desul::Impl::lock_address_cuda((void*)dest_val, scope)) {
                desul::atomic_thread_fence(desul::MemoryOrderAcquire(), scope);
                return_val = op.apply(*dest, val);
                *dest      = return_val;
                desul::atomic_thread_fence(desul::MemoryOrderRelease(), scope);
                desul::Impl::unlock_address_cuda((void*)dest_val, scope);
                done = 1;
              }
            }
            done_active = __ballot_sync(mask, done);
          }
          return return_val;
        }
      }

      template <typename Oper, typename DestPtrT, typename ValT, typename T>
      __device__
      typename Sacado::BaseExprType< Expr<T> >::type
      atomic_fetch_oper_device(const Oper& op, DestPtrT dest, ValT* dest_val,
                               const Expr<T>& x)
      {
        typedef typename Sacado::BaseExprType< Expr<T> >::type return_type;
        const typename Expr<T>::derived_type& val = x.derived();

        auto scope = desul::MemoryScopeDevice();

        if (atomics_use_team()) {
          int go = 1;
          while (go) {
            if (threadIdx.x == 0)
              go = !desul::Impl::lock_address_cuda((void*)dest_val, scope);
            go = Kokkos::shfl(go, 0, blockDim.x);
          }
          desul::atomic_thread_fence(desul::MemoryOrderAcquire(), scope);
          return_type return_val = *dest;
          *dest                  = op.apply(return_val, val);
          desul::atomic_thread_fence(desul::MemoryOrderRelease(), scope);
          if (threadIdx.x == 0)
            desul::Impl::unlock_address_cuda((void*)dest_val, scope);
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
              if (desul::Impl::lock_address_cuda((void*)dest_val, scope)) {
                desul::atomic_thread_fence(desul::MemoryOrderAcquire(), scope);
                return_val = *dest;
                *dest      = op.apply(return_val, val);
                desul::atomic_thread_fence(desul::MemoryOrderRelease(), scope);
                desul::Impl::unlock_address_cuda((void*)dest_val, scope);
                done = 1;
              }
            }
            done_active = __ballot_sync(mask, done);
          }
          return return_val;
        }
      }

#elif defined(KOKKOS_ENABLE_HIP)

      // Our implementation of Kokkos::atomic_oper_fetch() and
      // Kokkos::atomic_fetch_oper() for Sacado types on device
      template <typename Oper, typename DestPtrT, typename ValT, typename T>
      __device__
      typename Sacado::BaseExprType< Expr<T> >::type
      atomic_oper_fetch_device(const Oper& op, DestPtrT dest, ValT* dest_val,
                               const Expr<T>& x)
      {
        typedef typename Sacado::BaseExprType< Expr<T> >::type return_type;
        const typename Expr<T>::derived_type& val = x.derived();

        auto scope = desul::MemoryScopeDevice();

        if (atomics_use_team()) {
          int go = 1;
          while (go) {
            if (threadIdx.x == 0)
              go = !desul::Impl::lock_address_hip((void*)dest_val, scope);
            go = Kokkos::shfl(go, 0, blockDim.x);
          }
          desul::atomic_thread_fence(desul::MemoryOrderAcquire(), scope);
          return_type return_val = op.apply(*dest, val);
          *dest                  = return_val;
          desul::atomic_thread_fence(desul::MemoryOrderRelease(), scope);
          if (threadIdx.x == 0)
            desul::Impl::unlock_address_hip((void*)dest_val, scope);
          return return_val;
        }
        else {
          return_type return_val;
          int done                 = 0;
          unsigned int active      = __ballot(1);
          unsigned int done_active = 0;
          while (active != done_active) {
            if (!done) {
              if (desul::Impl::lock_address_hip((void*)dest_val, scope)) {
                return_val = op.apply(*dest, val);
                *dest      = return_val;
                desul::Impl::unlock_address_hip((void*)dest_val, scope);
                done = 1;
              }
            }
            done_active = __ballot(done);
          }
          return return_val;
        }
      }

      template <typename Oper, typename DestPtrT, typename ValT, typename T>
      __device__
      typename Sacado::BaseExprType< Expr<T> >::type
      atomic_fetch_oper_device(const Oper& op, DestPtrT dest, ValT* dest_val,
                               const Expr<T>& x)
      {
        typedef typename Sacado::BaseExprType< Expr<T> >::type return_type;
        const typename Expr<T>::derived_type& val = x.derived();

        auto scope = desul::MemoryScopeDevice();

        if (atomics_use_team()) {
          int go = 1;
          while (go) {
            if (threadIdx.x == 0)
              go = !desul::Impl::lock_address_hip((void*)dest_val, scope);
            go = Kokkos::shfl(go, 0, blockDim.x);
          }
          desul::atomic_thread_fence(desul::MemoryOrderAcquire(), scope);
          return_type return_val = *dest;
          *dest                  = op.apply(return_val, val);
          desul::atomic_thread_fence(desul::MemoryOrderRelease(), scope);
          if (threadIdx.x == 0)
            desul::Impl::unlock_address_hip((void*)dest_val, scope);
          return return_val;
        }
        else {
          return_type return_val;
          int done                 = 0;
          unsigned int active      = __ballot(1);
          unsigned int done_active = 0;
          while (active != done_active) {
            if (!done) {
              if (desul::Impl::lock_address_hip((void*)dest_val, scope)) {
                return_val = *dest;
                *dest      = op.apply(return_val, val);
                desul::Impl::unlock_address_hip((void*)dest_val, scope);
                done = 1;
              }
            }
            done_active = __ballot(done);
          }
          return return_val;
        }
      }

#elif defined(KOKKOS_ENABLE_SYCL)

      // Our implementation of Kokkos::atomic_oper_fetch() and
      // Kokkos::atomic_fetch_oper() for Sacado types on device
      template <typename Oper, typename DestPtrT, typename ValT, typename T>
      typename Sacado::BaseExprType< Expr<T> >::type
      atomic_oper_fetch_device(const Oper& op, DestPtrT dest, ValT* dest_val,
                               const Expr<T>& x)
      {
        Kokkos::abort("Not implemented!");
        return {};
      }

      template <typename Oper, typename DestPtrT, typename ValT, typename T>
      typename Sacado::BaseExprType< Expr<T> >::type
      atomic_fetch_oper_device(const Oper& op, DestPtrT dest, ValT* dest_val,
                               const Expr<T>& x)
      {
        Kokkos::abort("Not implemented!");
        return {};
      }
#endif

      // Overloads of Kokkos::atomic_oper_fetch/Kokkos::atomic_fetch_oper
      // for Sacado types
      template <typename Oper, typename S>
      SACADO_INLINE_FUNCTION GeneralFad<S>
      atomic_oper_fetch(const Oper& op, GeneralFad<S>* dest,
                        const GeneralFad<S>& val)
      {
        KOKKOS_IF_ON_HOST(return Impl::atomic_oper_fetch_host(op, dest, &(dest->val()), val);)
        KOKKOS_IF_ON_DEVICE(return Impl::atomic_oper_fetch_device(op, dest, &(dest->val()), val);)
      }
      template <typename Oper, typename ValT, unsigned sl, unsigned ss,
                typename U, typename T>
      SACADO_INLINE_FUNCTION U
      atomic_oper_fetch(const Oper& op, ViewFadPtr<ValT,sl,ss,U> dest,
                        const Expr<T>& val)
      {
        KOKKOS_IF_ON_HOST(return Impl::atomic_oper_fetch_host(op, dest, &dest.val(), val);)
        KOKKOS_IF_ON_DEVICE(return Impl::atomic_oper_fetch_device(op, dest, &dest.val(), val);)
      }

      template <typename Oper, typename S>
      SACADO_INLINE_FUNCTION GeneralFad<S>
      atomic_fetch_oper(const Oper& op, GeneralFad<S>* dest,
                        const GeneralFad<S>& val)
      {
        KOKKOS_IF_ON_HOST(return Impl::atomic_fetch_oper_host(op, dest, &(dest->val()), val);)
        KOKKOS_IF_ON_DEVICE(return Impl::atomic_fetch_oper_device(op, dest, &(dest->val()), val);)
      }
      template <typename Oper, typename ValT, unsigned sl, unsigned ss,
                typename U, typename T>
      SACADO_INLINE_FUNCTION U
      atomic_fetch_oper(const Oper& op, ViewFadPtr<ValT,sl,ss,U> dest,
                        const Expr<T>& val)
      {
        KOKKOS_IF_ON_HOST(return Impl::atomic_fetch_oper_host(op, dest, &dest.val(), val);)
        KOKKOS_IF_ON_DEVICE(return Impl::atomic_fetch_oper_device(op, dest, &dest.val(), val);)
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

#endif // HAVE_SACADO_KOKKOS
#endif // SACADO_FAD_EXP_VIEWFAD_HPP
