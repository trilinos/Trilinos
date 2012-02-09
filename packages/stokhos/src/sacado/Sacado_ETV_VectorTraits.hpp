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

#ifndef SACADO_ETV_VECTORTRAITS_HPP
#define SACADO_ETV_VECTORTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace ETV {
    template <typename T, typename S> class Vector;
  }
}

namespace Sacado {

  //! Specialization of %Promote to Taylor types
  template <typename T, typename S>
  class Promote< ETV::Vector<T,S>, ETV::Vector<T,S> > {
  public:

    typedef ETV::Vector<T,S> type;
  };

  //! Specialization of %Promote to Vector types
  template <typename L, typename R, typename S>
  class Promote< ETV::Vector<L,S>, R > {
  public:

    typedef typename ValueType< ETV::Vector<L,S> >::type value_type_l;
    typedef typename ValueType<R>::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef ETV::Vector<value_type,S> type;
  };

  //! Specialization of %Promote to Vector types
  template <typename L, typename R, typename S>
  class Promote< L, ETV::Vector<R,S> > {
  public:

    typedef typename ValueType<L>::type value_type_l;
    typedef typename ValueType< ETV::Vector<R,S> >::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef ETV::Vector<value_type,S> type;
  };

  //! Specialization of %ScalarType to Vector types
  template <typename T, typename S>
  struct ScalarType< ETV::Vector<T,S> > {
    typedef typename ScalarType<typename ETV::Vector<T,S>::value_type>::type type;
  };

  //! Specialization of %ValueType to Vector types
  template <typename T, typename S>
  struct ValueType< ETV::Vector<T,S> > {
    typedef typename ETV::Vector<T,S>::value_type type;
  };

  //! Specialization of %IsADType to Vector types
  template <typename T, typename S>
  struct IsADType< ETV::Vector<T,S> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to Vector types
  template <typename T, typename S>
  struct IsScalarType< ETV::Vector<T,S> > {
    static const bool value = false;
  };

  //! Specialization of %Value to Vector types
  template <typename T, typename S>
  struct Value< ETV::Vector<T,S> > {
    typedef typename ValueType< ETV::Vector<T,S> >::type value_type;
    static const value_type& eval(const ETV::Vector<T,S>& x) { 
      return x.val(); }
  };

  //! Specialization of %ScalarValue to Vector types
  template <typename T, typename S>
  struct ScalarValue< ETV::Vector<T,S> > {
    typedef typename ValueType< ETV::Vector<T,S> >::type value_type;
    typedef typename ScalarType< ETV::Vector<T,S> >::type scalar_type;
    static const scalar_type& eval(const ETV::Vector<T,S>& x) { 
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to Vector types
  template <typename T, typename S>
  struct StringName< ETV::Vector<T,S> > {
    static std::string eval() { 
      return std::string("Sacado::ETV::Vector< ") + 
	StringName<T>::eval() + " >"; }
  };

  //! Specialization of IsEqual to Vector types
  template <typename T, typename S>
  struct IsEqual< ETV::Vector<T,S> > {
    static bool eval(const ETV::Vector<T,S>& x, 
		     const ETV::Vector<T,S>& y) {
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
  struct PromotionTraits< Sacado::ETV::Vector<T,S>, 
			  Sacado::ETV::Vector<T,S> > {
    typedef typename Sacado::Promote< Sacado::ETV::Vector<T,S>,
				      Sacado::ETV::Vector<T,S> >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename T, typename S, typename R>
  struct PromotionTraits< Sacado::ETV::Vector<T,S>, R > {
    typedef typename Sacado::Promote< Sacado::ETV::Vector<T,S>, R >::type 
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename L, typename T, typename S>
  struct PromotionTraits< L, Sacado::ETV::Vector<T,S> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::ETV::Vector<T,S> >::type 
    promote;
  };

  //! Specializtion of Teuchos::ScalarTraits
  template <typename T, typename S>
  struct ScalarTraits< Sacado::ETV::Vector<T,S> > : 
    public Sacado::ETV::ScalarTraitsImp< Sacado::ETV::Vector<T,S> > {};


  //! Specialization of %Teuchos::SerializationTraits
  template <typename Ordinal, typename T, typename S>
  struct SerializationTraits<Ordinal, Sacado::ETV::Vector<T,S> > :
    public Sacado::ETV::SerializationTraitsImp< Ordinal, 
						Sacado::ETV::Vector<T,S> > {};

  //! Specialization of %Teuchos::ValueTypeSerializer
  template <typename Ordinal, typename T, typename S>
  struct ValueTypeSerializer<Ordinal, Sacado::ETV::Vector<T,S> > :
    public Sacado::ETV::SerializerImp< Ordinal,
				       Sacado::ETV::Vector<T,S>,
				       ValueTypeSerializer<Ordinal,T> >
  {
    typedef Sacado::ETV::Vector<T,S> VecType;
    typedef ValueTypeSerializer<Ordinal,T> ValueSerializer;
    typedef Sacado::ETV::SerializerImp< Ordinal,VecType,ValueSerializer> Base;
    ValueTypeSerializer(const Teuchos::RCP<const ValueSerializer>& vs,
			Ordinal sz = 0) :
      Base(vs, sz) {}
  };
  
}
#endif // HAVE_SACADO_TEUCHOS

#endif // SACADO_ETV_VECTORTRAITS_HPP
