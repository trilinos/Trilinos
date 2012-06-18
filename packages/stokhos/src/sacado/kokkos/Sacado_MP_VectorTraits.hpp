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

#ifndef SACADO_MP_VECTOR_TRAITS_HPP
#define SACADO_MP_VECTOR_TRAITS_HPP

#include "Sacado_Traits.hpp"
#include "Sacado_mpl_apply.hpp"

//
// Currently only the type-style traits classes will work on the device
// since the function-style traits don't have the needed specializations
// for the device.  This is in particular true for the scalar traits.
//

// Forward declarations
namespace Sacado {
  namespace MP {
    template <typename S, typename N> class Vector;
  }
}

namespace Sacado {

  //! Specialization of %Promote to Taylor types
  template <typename S, typename N>
  class Promote< MP::Vector<S,N>, MP::Vector<S,N> > {
  public:

    typedef MP::Vector<S,N> type;
  };

  //! Specialization of %Promote to Vector types
  template <typename R, typename S, typename N>
  class Promote< MP::Vector<S,N>, R > {
  public:

    typedef typename ValueType< MP::Vector<S,N> >::type value_type_l;
    typedef typename ValueType<R>::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;
    typedef typename S::ordinal_type ordinal_type;
    typedef typename Sacado::mpl::apply<S,ordinal_type,value_type>::type storage_type;

    typedef MP::Vector<storage_type,N> type;
  };

  //! Specialization of %Promote to Vector types
  template <typename L, typename S, typename N>
  class Promote< L, MP::Vector<S,N> > {
  public:

    typedef typename ValueType<L>::type value_type_l;
    typedef typename ValueType< MP::Vector<S,N> >::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;
    typedef typename S::ordinal_type ordinal_type;
    typedef typename Sacado::mpl::apply<S,ordinal_type,value_type>::type storage_type;

    typedef MP::Vector<storage_type,N> type;
  };

  //! Specialization of %ScalarType to Vector types
  template <typename S, typename N>
  struct ScalarType< MP::Vector<S,N> > {
    typedef typename ScalarType<typename MP::Vector<S,N>::value_type>::type type;
  };

  //! Specialization of %ValueType to Vector types
  template <typename S, typename N>
  struct ValueType< MP::Vector<S,N> > {
    typedef typename MP::Vector<S,N>::value_type type;
  };

  //! Specialization of %IsADType to Vector types
  template <typename S, typename N>
  struct IsADType< MP::Vector<S,N> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to Vector types
  template <typename S, typename N>
  struct IsScalarType< MP::Vector<S,N> > {
    static const bool value = false;
  };

  //! Specialization of %Value to Vector types
  template <typename S, typename N>
  struct Value< MP::Vector<S,N> > {
    typedef typename ValueType< MP::Vector<S,N> >::type value_type;
    static const value_type& eval(const MP::Vector<S,N>& x) { 
      return x.val(); }
  };

  //! Specialization of %ScalarValue to Vector types
  template <typename S, typename N>
  struct ScalarValue< MP::Vector<S,N> > {
    typedef typename ValueType< MP::Vector<S,N> >::type value_type;
    typedef typename ScalarType< MP::Vector<S,N> >::type scalar_type;
    static const scalar_type& eval(const MP::Vector<S,N>& x) { 
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to Vector types
  template <typename S, typename N>
  struct StringName< MP::Vector<S,N> > {
    static std::string eval() { 
      return std::string("Sacado::MP::Vector< ") + 
	StringName<S>::eval() + " >"; }
  };

  //! Specialization of IsEqual to Vector types
  template <typename S, typename N>
  struct IsEqual< MP::Vector<S,N> > {
    static bool eval(const MP::Vector<S,N>& x, 
		     const MP::Vector<S,N>& y) {
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
  template <typename S, typename N>
  struct PromotionTraits< Sacado::MP::Vector<S,N>, 
			  Sacado::MP::Vector<S,N> > {
    typedef typename Sacado::Promote< Sacado::MP::Vector<S,N>,
				      Sacado::MP::Vector<S,N> >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename S, typename N, typename R>
  struct PromotionTraits< Sacado::MP::Vector<S,N>, R > {
    typedef typename Sacado::Promote< Sacado::MP::Vector<S,N>, R >::type 
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename L, typename S, typename N>
  struct PromotionTraits< L, Sacado::MP::Vector<S,N> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::MP::Vector<S,N> >::type 
    promote;
  };

  //! Specializtion of Teuchos::ScalarTraits
  template <typename S, typename N>
  struct ScalarTraits< Sacado::MP::Vector<S,N> > : 
    public Sacado::ETV::ScalarTraitsImp< Sacado::MP::Vector<S,N> > {};


  //! Specialization of %Teuchos::SerializationTraits
  template <typename Ordinal, typename S, typename N>
  struct SerializationTraits<Ordinal, Sacado::MP::Vector<S,N> > :
    public Sacado::ETV::SerializationTraitsImp< Ordinal, 
						Sacado::MP::Vector<S,N> > {};

  //! Specialization of %Teuchos::ValueTypeSerializer
  template <typename Ordinal, typename S, typename N>
  struct ValueTypeSerializer<Ordinal, Sacado::MP::Vector<S,N> > :
    public Sacado::ETV::SerializerImp< Ordinal,
				       Sacado::MP::Vector<S,N>,
				       ValueTypeSerializer<Ordinal,typename Sacado::MP::Vector<S,N>::value_type> >
  {
    typedef Sacado::MP::Vector<S,N> VecType;
    typedef typename VecType::value_type value_type;
    typedef ValueTypeSerializer<Ordinal,value_type> ValueSerializer;
    typedef Sacado::ETV::SerializerImp< Ordinal,VecType,ValueSerializer> Base;
    ValueTypeSerializer(const Teuchos::RCP<const ValueSerializer>& vs,
			Ordinal sz = 0) :
      Base(vs, sz) {}
  };
  
}
#endif // HAVE_SACADO_TEUCHOS

#endif // SACADO_MP_VECTORTRAITS_HPP
