// $Id$ 
// $Source$ 
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

#ifndef SACADO_FAD_DMFADTRAITS_HPP
#define SACADO_FAD_DMFADTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace Fad {
    template <typename T> class DMFad;
  }
}

namespace Sacado {

  //! Specialization of %Promote to DMFad types
  SACADO_FAD_PROMOTE_SPEC( Fad, DMFad )

  //! Specialization of %ScalarType to DMFad types
  template <typename ValueT>
  struct ScalarType< Fad::DMFad<ValueT> > {
    typedef typename Fad::DMFad<ValueT>::ScalarT type;
  };

  //! Specialization of %ValueType to DMFad types
  template <typename ValueT>
  struct ValueType< Fad::DMFad<ValueT> > {
    typedef ValueT type;
  };

  //! Specialization of %IsADType to DMFad types
  template <typename ValueT>
  struct IsADType< Fad::DMFad<ValueT> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to DMFad types
  template <typename ValueT>
  struct IsScalarType< Fad::DMFad<ValueT> > {
    static const bool value = false;
  };

  //! Specialization of %Value to DMFad types
  template <typename ValueT>
  struct Value< Fad::DMFad<ValueT> > {
    typedef typename ValueType< Fad::DMFad<ValueT> >::type value_type;
    static const value_type& eval(const Fad::DMFad<ValueT>& x) { 
      return x.val(); }
  };

  //! Specialization of %ScalarValue to DMFad types
  template <typename ValueT>
  struct ScalarValue< Fad::DMFad<ValueT> > {
    typedef typename ValueType< Fad::DMFad<ValueT> >::type value_type;
    typedef typename ScalarType< Fad::DMFad<ValueT> >::type scalar_type;
    static const scalar_type& eval(const Fad::DMFad<ValueT>& x) { 
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to DMFad types
  template <typename ValueT>
  struct StringName< Fad::DMFad<ValueT> > {
    static std::string eval() { 
      return std::string("Sacado::Fad::DMFad< ") + 
	StringName<ValueT>::eval() + " >"; }
  };

  //! Specialization of %IsEqual to DMFad types
  template <typename ValueT>
  struct IsEqual< Fad::DMFad<ValueT> > {
    static bool eval(const Fad::DMFad<ValueT>& x, const Fad::DMFad<ValueT>& y) {
      return x.isEqualTo(y);
    }
  };

  //! Specialization of %IsStaticallySized to DMFad types
  template <typename ValueT>
  struct IsStaticallySized< Fad::DMFad<ValueT> > {
    static const bool value = false;
  };

  template <typename T>
  struct IsFad< Fad::DMFad<T> > {
    static const bool value = true;
  };

} // namespace Sacado

// Define Teuchos traits classes
#ifdef HAVE_SACADO_TEUCHOS
#include "Teuchos_PromotionTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Sacado_Fad_ScalarTraitsImp.hpp"

namespace Teuchos {

  //! Specialization of %Teuchos::PromotionTraits to DMFad types
  template <typename ValueT>
  struct PromotionTraits< Sacado::Fad::DMFad<ValueT>, 
			  Sacado::Fad::DMFad<ValueT> > {
    typedef typename Sacado::Promote< Sacado::Fad::DMFad<ValueT>,
				      Sacado::Fad::DMFad<ValueT> >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DMFad types
  template <typename ValueT, typename R>
  struct PromotionTraits< Sacado::Fad::DMFad<ValueT>, R > {
    typedef typename Sacado::Promote< Sacado::Fad::DMFad<ValueT>, R >::type 
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DMFad types
  template <typename L, typename ValueT>
  struct PromotionTraits< L, Sacado::Fad::DMFad<ValueT> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::Fad::DMFad<ValueT> >::type 
    promote;
  };

  //! Specializtion of %Teuchos::ScalarTraits
  template <typename ValueT>
  struct ScalarTraits< Sacado::Fad::DMFad<ValueT> > :
    public Sacado::Fad::ScalarTraitsImp< Sacado::Fad::DMFad<ValueT> >
  {};

  //! Specialization of %Teuchos::SerializationTraits
  template <typename Ordinal, typename ValueT>
  struct SerializationTraits<Ordinal, Sacado::Fad::DMFad<ValueT> > :
    public Sacado::Fad::SerializationTraitsImp< Ordinal, 
						Sacado::Fad::DMFad<ValueT> > 
  {};

  //! Specialization of %Teuchos::ValueTypeSerializer
  template <typename Ordinal, typename ValueT>
  struct ValueTypeSerializer<Ordinal, Sacado::Fad::DMFad<ValueT> > :
    public Sacado::Fad::SerializerImp< Ordinal, 
				       Sacado::Fad::DMFad<ValueT>,
				       ValueTypeSerializer<Ordinal,ValueT> > 
  {
    typedef Sacado::Fad::DMFad<ValueT> FadType;
    typedef ValueTypeSerializer<Ordinal,ValueT> ValueSerializer;
    typedef Sacado::Fad::SerializerImp< Ordinal,FadType,ValueSerializer> Base;
    ValueTypeSerializer(const Teuchos::RCP<const ValueSerializer>& vs,
			Ordinal sz = 0) :
      Base(vs, sz) {}
  };
}
#endif // HAVE_SACADO_TEUCHOS

#endif // SACADO_FAD_DMFADTRAITS_HPP
