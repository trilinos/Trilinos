// @HEADER
// *****************************************************************************
//                           Sacado Package
//
// Copyright 2006 NTESS and the Sacado contributors.
// SPDX-License-Identifier: LGPL-2.1-or-later
//
// ***********************************************************************
//
// The forward-mode AD classes in Sacado are a derivative work of the
// expression template classes in the Fad package by Nicolas Di Cesare.
// The following banner is included in the original Fad source code:
//
// ************ DO NOT REMOVE THIS BANNER ****************
//
//  Nicolas Di Cesare <Nicolas.Dicesare@ann.jussieu.fr>
//  http://www.ann.jussieu.fr/~dicesare
//
//            CEMRACS 98 : C++ courses,
//         templates : new C++ techniques
//            for scientific computing
//
//********************************************************
//
//  NumericalTraits class to illustrate TRAITS
//
//********************************************************
// @HEADER

#ifndef SACADO_CACHEFAD_VIEWFADTRAITS_HPP
#define SACADO_CACHEFAD_VIEWFADTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace CacheFad {
    template <typename T,unsigned,unsigned,typename> class ViewFad;
  }
}

namespace Sacado {

  //! Specialization of %Promote to ViewFad types
  SACADO_VFAD_PROMOTE_SPEC( CacheFad )

  //! Specialization of %ScalarType to ViewFad types
  template <typename ValueT, unsigned Size, unsigned Stride, typename Base>
  struct ScalarType< CacheFad::ViewFad<ValueT,Size,Stride,Base> > {
    typedef typename CacheFad::ViewFad<ValueT,Size,Stride,Base>::ScalarT type;
  };

  //! Specialization of %ValueType to ViewFad types
  template <typename ValueT, unsigned Size, unsigned Stride, typename Base>
  struct ValueType< CacheFad::ViewFad<ValueT,Size,Stride,Base> > {
    typedef ValueT type;
  };

  //! Specialization of %IsADType to ViewFad types
  template <typename ValueT, unsigned Size, unsigned Stride, typename Base>
  struct IsADType< CacheFad::ViewFad<ValueT,Size,Stride,Base> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to ViewFad types
  template <typename ValueT, unsigned Size, unsigned Stride, typename Base>
  struct IsScalarType< CacheFad::ViewFad<ValueT,Size,Stride,Base> > {
    static const bool value = false;
  };

  //! Specialization of %Value to ViewFad types
  template <typename ValueT, unsigned Size, unsigned Stride, typename Base>
  struct Value< CacheFad::ViewFad<ValueT,Size,Stride,Base> > {
    typedef typename ValueType< CacheFad::ViewFad<ValueT,Size,Stride,Base> >::type value_type;
    SACADO_INLINE_FUNCTION
    static const value_type& eval(const CacheFad::ViewFad<ValueT,Size,Stride,Base>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to ViewFad types
  template <typename ValueT, unsigned Size, unsigned Stride, typename Base>
  struct ScalarValue< CacheFad::ViewFad<ValueT,Size,Stride,Base> > {
    typedef typename ValueType< CacheFad::ViewFad<ValueT,Size,Stride,Base> >::type value_type;
    typedef typename ScalarType< CacheFad::ViewFad<ValueT,Size,Stride,Base> >::type scalar_type;
    SACADO_INLINE_FUNCTION
    static const scalar_type& eval(const CacheFad::ViewFad<ValueT,Size,Stride,Base>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to ViewFad types
  template <typename ValueT, unsigned Size, unsigned Stride, typename Base>
  struct StringName< CacheFad::ViewFad<ValueT,Size,Stride,Base> > {
    static std::string eval() {
      return std::string("Sacado::CacheFad::ViewFad< ") +
        StringName<ValueT>::eval() + " >"; }
  };

  //! Specialization of %IsEqual to ViewFad types
  template <typename ValueT, unsigned Size, unsigned Stride, typename Base>
  struct IsEqual< CacheFad::ViewFad<ValueT,Size,Stride,Base> > {
    SACADO_INLINE_FUNCTION
    static bool eval(const CacheFad::ViewFad<ValueT,Size,Stride,Base>& x,
                     const CacheFad::ViewFad<ValueT,Size,Stride,Base>& y) {
      return x.isEqualTo(y);
    }
  };

  //! Specialization of %IsStaticallySized to ViewFad types
  template <typename ValueT, unsigned Size, unsigned Stride, typename Base>
  struct IsStaticallySized< CacheFad::ViewFad<ValueT,Size,Stride,Base> > {
    static const bool value = false;
  };

} // namespace Sacado

// ViewFad is not a proper scalar type, so we don't define any of the
// Teuchos traits classes

#endif // SACADO_CACHEFAD_VIEWFADTRAITS_HPP
