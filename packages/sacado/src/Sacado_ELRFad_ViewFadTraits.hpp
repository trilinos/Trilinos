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

#ifndef SACADO_ELRFAD_VIEWFADTRAITS_HPP
#define SACADO_ELRFAD_VIEWFADTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace ELRFad {
    template <typename T,unsigned,unsigned,typename> class ViewFad;
  }
}

namespace Sacado {

  //! Specialization of %Promote to ViewFad types
  SACADO_VFAD_PROMOTE_SPEC( ELRFad )

  //! Specialization of %ScalarType to ViewFad types
  template <typename ValueT, unsigned Size, unsigned Stride, typename Base>
  struct ScalarType< ELRFad::ViewFad<ValueT,Size,Stride,Base> > {
    typedef typename ELRFad::ViewFad<ValueT,Size,Stride,Base>::ScalarT type;
  };

  //! Specialization of %ValueType to ViewFad types
  template <typename ValueT, unsigned Size, unsigned Stride, typename Base>
  struct ValueType< ELRFad::ViewFad<ValueT,Size,Stride,Base> > {
    typedef ValueT type;
  };

  //! Specialization of %IsADType to ViewFad types
  template <typename ValueT, unsigned Size, unsigned Stride, typename Base>
  struct IsADType< ELRFad::ViewFad<ValueT,Size,Stride,Base> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to ViewFad types
  template <typename ValueT, unsigned Size, unsigned Stride, typename Base>
  struct IsScalarType< ELRFad::ViewFad<ValueT,Size,Stride,Base> > {
    static const bool value = false;
  };

  //! Specialization of %Value to ViewFad types
  template <typename ValueT, unsigned Size, unsigned Stride, typename Base>
  struct Value< ELRFad::ViewFad<ValueT,Size,Stride,Base> > {
    typedef typename ValueType< ELRFad::ViewFad<ValueT,Size,Stride,Base> >::type value_type;
    KOKKOS_INLINE_FUNCTION
    static const value_type& eval(const ELRFad::ViewFad<ValueT,Size,Stride,Base>& x) {
      return x.val(); }
  };

  //! Specialization of %ScalarValue to ViewFad types
  template <typename ValueT, unsigned Size, unsigned Stride, typename Base>
  struct ScalarValue< ELRFad::ViewFad<ValueT,Size,Stride,Base> > {
    typedef typename ValueType< ELRFad::ViewFad<ValueT,Size,Stride,Base> >::type value_type;
    typedef typename ScalarType< ELRFad::ViewFad<ValueT,Size,Stride,Base> >::type scalar_type;
    KOKKOS_INLINE_FUNCTION
    static const scalar_type& eval(const ELRFad::ViewFad<ValueT,Size,Stride,Base>& x) {
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to ViewFad types
  template <typename ValueT, unsigned Size, unsigned Stride, typename Base>
  struct StringName< ELRFad::ViewFad<ValueT,Size,Stride,Base> > {
    static std::string eval() {
      return std::string("Sacado::ELRFad::ViewFad< ") +
        StringName<ValueT>::eval() + " >"; }
  };

  //! Specialization of %IsEqual to ViewFad types
  template <typename ValueT, unsigned Size, unsigned Stride, typename Base>
  struct IsEqual< ELRFad::ViewFad<ValueT,Size,Stride,Base> > {
    KOKKOS_INLINE_FUNCTION
    static bool eval(const ELRFad::ViewFad<ValueT,Size,Stride,Base>& x,
                     const ELRFad::ViewFad<ValueT,Size,Stride,Base>& y) {
      return x.isEqualTo(y);
    }
  };

  //! Specialization of %IsStaticallySized to ViewFad types
  template <typename ValueT, unsigned Size, unsigned Stride, typename Base>
  struct IsStaticallySized< ELRFad::ViewFad<ValueT,Size,Stride,Base> > {
    static const bool value = false;
  };

} // namespace Sacado

// ViewFad is not a proper scalar type, so we don't define any of the
// Teuchos traits classes

#endif // SACADO_ELRFAD_VIEWFADTRAITS_HPP
