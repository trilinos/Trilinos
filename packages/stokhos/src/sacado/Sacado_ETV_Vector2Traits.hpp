// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_ETV_VECTOR2TRAITS_HPP
#define SACADO_ETV_VECTOR2TRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace ETV {
    template <typename T, typename S> class Vector2;
  }
}

namespace Sacado {

  //! Specialization of %Promote to Taylor types
  template <typename T, typename S>
  class Promote< ETV::Vector2<T,S>, ETV::Vector2<T,S> > {
  public:

    typedef ETV::Vector2<T,S> type;
  };

  //! Specialization of %Promote to Vector2 types
  template <typename L, typename R, typename S>
  class Promote< ETV::Vector2<L,S>, R > {
  public:

    typedef typename ValueType< ETV::Vector2<L,S> >::type value_type_l;
    typedef typename ValueType<R>::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef ETV::Vector2<value_type,S> type;
  };

  //! Specialization of %Promote to Vector2 types
  template <typename L, typename R, typename S>
  class Promote< L, ETV::Vector2<R,S> > {
  public:

    typedef typename ValueType<L>::type value_type_l;
    typedef typename ValueType< ETV::Vector2<R,S> >::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef ETV::Vector2<value_type,S> type;
  };

  //! Specialization of %ScalarType to Vector2 types
  template <typename T, typename S>
  struct ScalarType< ETV::Vector2<T,S> > {
    typedef typename ScalarType<typename ETV::Vector2<T,S>::value_type>::type type;
  };

  //! Specialization of %ValueType to Vector2 types
  template <typename T, typename S>
  struct ValueType< ETV::Vector2<T,S> > {
    typedef typename ETV::Vector2<T,S>::value_type type;
  };

  //! Specialization of %IsADType to Vector2 types
  template <typename T, typename S>
  struct IsADType< ETV::Vector2<T,S> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to Vector2 types
  template <typename T, typename S>
  struct IsScalarType< ETV::Vector2<T,S> > {
    static const bool value = false;
  };

  //! Specialization of %Value to Vector2 types
  template <typename T, typename S>
  struct Value< ETV::Vector2<T,S> > {
    typedef typename ValueType< ETV::Vector2<T,S> >::type value_type;
    static const value_type& eval(const ETV::Vector2<T,S>& x) { 
      return x.val(); }
  };

  //! Specialization of %ScalarValue to Vector2 types
  template <typename T, typename S>
  struct ScalarValue< ETV::Vector2<T,S> > {
    typedef typename ValueType< ETV::Vector2<T,S> >::type value_type;
    typedef typename ScalarType< ETV::Vector2<T,S> >::type scalar_type;
    static const scalar_type& eval(const ETV::Vector2<T,S>& x) { 
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to Vector2 types
  template <typename T, typename S>
  struct StringName< ETV::Vector2<T,S> > {
    static std::string eval() { 
      return std::string("Sacado::ETV::Vector2< ") + 
	StringName<T>::eval() + " >"; }
  };

  //! Specialization of IsEqual to Vector2 types
  template <typename T, typename S>
  struct IsEqual< ETV::Vector2<T,S> > {
    static bool eval(const ETV::Vector2<T,S>& x, 
		     const ETV::Vector2<T,S>& y) {
      return x.isEqualTo(y);
    }
  };

} // namespace Sacado

// Define Teuchos traits classes
#ifdef HAVE_SACADO_TEUCHOS
#include "Teuchos_PromotionTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Sacado_ETV_ScalarTraitsImp.hpp"
#include "Teuchos_SerializationTraits.hpp"

namespace Teuchos {

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename T, typename S>
  struct PromotionTraits< Sacado::ETV::Vector2<T,S>, 
			  Sacado::ETV::Vector2<T,S> > {
    typedef typename Sacado::Promote< Sacado::ETV::Vector2<T,S>,
				      Sacado::ETV::Vector2<T,S> >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename T, typename S, typename R>
  struct PromotionTraits< Sacado::ETV::Vector2<T,S>, R > {
    typedef typename Sacado::Promote< Sacado::ETV::Vector2<T,S>, R >::type 
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename L, typename T, typename S>
  struct PromotionTraits< L, Sacado::ETV::Vector2<T,S> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::ETV::Vector2<T,S> >::type 
    promote;
  };

  //! Specializtion of Teuchos::ScalarTraits
  template <typename T, typename S>
  struct ScalarTraits< Sacado::ETV::Vector2<T,S> > : 
    public Sacado::ETV::ScalarTraitsImp< Sacado::ETV::Vector2<T,S> > {};


  //! Specialization of %Teuchos::SerializationTraits
  template <typename Ordinal, typename T, typename S>
  struct SerializationTraits<Ordinal, Sacado::ETV::Vector2<T,S> > :
    public Sacado::ETV::SerializationTraitsImp< Ordinal, 
						Sacado::ETV::Vector2<T,S> > {};

  //! Specialization of %Teuchos::ValueTypeSerializer
  template <typename Ordinal, typename T, typename S>
  struct ValueTypeSerializer<Ordinal, Sacado::ETV::Vector2<T,S> > :
    public Sacado::ETV::SerializerImp< Ordinal,
				       Sacado::ETV::Vector2<T,S>,
				       ValueTypeSerializer<Ordinal,T> >
  {
    typedef Sacado::ETV::Vector2<T,S> VecType;
    typedef ValueTypeSerializer<Ordinal,T> ValueSerializer;
    typedef Sacado::ETV::SerializerImp< Ordinal,VecType,ValueSerializer> Base;
    ValueTypeSerializer(const Teuchos::RCP<const ValueSerializer>& vs,
			Ordinal sz = 0) :
      Base(vs, sz) {}
  };
  
}
#endif // HAVE_SACADO_TEUCHOS

#endif // SACADO_ETV_VECTORTRAITS_HPP
