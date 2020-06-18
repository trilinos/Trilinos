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

#ifndef SACADO_FAD_EXP_VIEWFAD_HPP
#define SACADO_FAD_EXP_VIEWFAD_HPP

#include "Sacado_Fad_Exp_GeneralFad.hpp"
#include "Sacado_Fad_Exp_ViewStorage.hpp"

#if defined(HAVE_SACADO_KOKKOSCORE)
#include "Kokkos_Atomic.hpp"
#include "impl/Kokkos_Error.hpp"
#endif

namespace Sacado {

  namespace Fad {
  namespace Exp {

    template <typename T, unsigned static_length, unsigned static_stride, typename U>
    using  ViewFad = GeneralFad< ViewStorage<T,static_length,static_stride,U> >;

    // Class representing a pointer to ViewFad so that &ViewFad is supported
    template <typename T, unsigned sl, unsigned ss, typename U>
    class ViewFadPtr : public ViewFad<T,sl,ss,U> {
    public:

      // Storage type base class
      typedef ViewFad<T,sl,ss,U> view_fad_type;

      // Bring in constructors
      using view_fad_type::view_fad_type;

      // Add overload of dereference operator
      KOKKOS_INLINE_FUNCTION
      view_fad_type* operator->() { return this; }

      // Add overload of dereference operator
      KOKKOS_INLINE_FUNCTION
      view_fad_type& operator*() { *this; }
    };

#if defined(HAVE_SACADO_KOKKOSCORE)
    // Overload of Kokkos::atomic_add for ViewFad types.
    template <typename ValT, unsigned sl, unsigned ss, typename U, typename T>
    KOKKOS_INLINE_FUNCTION
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
#endif

  } // namespace Exp
  } // namespace Fad

  template <typename,unsigned,unsigned> struct ViewFadType;

  //! The View Fad type associated with this type
  template< class S, unsigned length, unsigned stride >
  struct ViewFadType< Fad::Exp::GeneralFad<S>, length, stride > {
    typedef Fad::Exp::ViewFad< typename S::value_type,length,stride,Fad::Exp::GeneralFad<S> > type;
  };

  //! The View Fad type associated with this type
  /*!
   * Do not include the const in the base expr type.
   */
  template< class S, unsigned length, unsigned stride >
  struct ViewFadType< const Fad::Exp::GeneralFad<S>, length, stride > {
    typedef Fad::Exp::ViewFad< const typename S::value_type,length,stride,Fad::Exp::GeneralFad<S> > type;
  };

  // Specialization of BaseExprType for ViewFad, to use the base fad type
  template <typename T, unsigned static_length, unsigned static_stride, typename U>
  struct BaseExprType< Fad::Exp::GeneralFad< Fad::Exp::ViewStorage<T,static_length,static_stride,U> > > {
    typedef U type;
  };

  //! Specialization of %ScalarType to ViewFad types
  /*!
   * This specialization overrides the one for GeneralFad to handle const
   * value types so that resulting scalar type is still const.
   */
  template <typename ValueT, unsigned Size, unsigned Stride, typename Base>
  struct ScalarType< Fad::Exp::ViewFad<ValueT,Size,Stride,Base> > {
    typedef typename ScalarType<ValueT>::type type;
  };

  /*!
   * This specialization overrides the one for GeneralFad to handle const
   * value types so that resulting value type is still const.
   */
  template <typename ValueT, unsigned Size, unsigned Stride, typename Base>
  struct ValueType< Fad::Exp::ViewFad<ValueT,Size,Stride,Base> > {
    typedef ValueT type;
  };

} // namespace Sacado

#endif // SACADO_FAD_EXP_VIEWFAD_HPP
