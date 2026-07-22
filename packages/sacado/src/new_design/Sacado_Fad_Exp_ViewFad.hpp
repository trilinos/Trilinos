// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
// *****************************************************************************
// @HEADER

#ifndef SACADO_FAD_EXP_VIEWFAD_HPP
#define SACADO_FAD_EXP_VIEWFAD_HPP

#include "Sacado_Fad_Exp_GeneralFad.hpp"
#include "Sacado_Fad_Exp_ViewStorage.hpp"

#if defined(HAVE_SACADO_KOKKOS)
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
      SACADO_INLINE_FUNCTION
      view_fad_type* operator->() { return this; }

      // Add overload of dereference operator
      SACADO_INLINE_FUNCTION
      view_fad_type& operator*() { return *this; }
    };

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
