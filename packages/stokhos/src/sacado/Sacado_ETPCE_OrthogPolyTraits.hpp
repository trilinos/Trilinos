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

#ifndef SACADO_ETPCE_ORTHOGPOLYTRAITS_HPP
#define SACADO_ETPCE_ORTHOGPOLYTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace ETPCE {
    template <typename T, typename S> class OrthogPoly;
  }
}

namespace Sacado {

  //! Specialization of %Promote to Taylor types
  template <typename T, typename S>
  class Promote< ETPCE::OrthogPoly<T,S>, ETPCE::OrthogPoly<T,S> > {
  public:

    typedef ETPCE::OrthogPoly<T,S> type;
  };

  //! Specialization of %Promote to OrthogPoly types
  template <typename L, typename R, typename S>
  class Promote< ETPCE::OrthogPoly<L,S>, R > {
  public:

    typedef typename ValueType< ETPCE::OrthogPoly<L,S> >::type value_type_l;
    typedef typename ValueType<R>::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef ETPCE::OrthogPoly<value_type,S> type;
  };

  //! Specialization of %Promote to OrthogPoly types
  template <typename L, typename R, typename S>
  class Promote< L, ETPCE::OrthogPoly<R,S> > {
  public:

    typedef typename ValueType<L>::type value_type_l;
    typedef typename ValueType< ETPCE::OrthogPoly<R,S> >::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef ETPCE::OrthogPoly<value_type,S> type;
  };

  //! Specialization of %ScalarType to OrthogPoly types
  template <typename T, typename S>
  struct ScalarType< ETPCE::OrthogPoly<T,S> > {
    typedef typename ScalarType<typename ETPCE::OrthogPoly<T,S>::value_type>::type type;
  };

  //! Specialization of %ValueType to OrthogPoly types
  template <typename T, typename S>
  struct ValueType< ETPCE::OrthogPoly<T,S> > {
    typedef typename ETPCE::OrthogPoly<T,S>::value_type type;
  };

  //! Specialization of %IsADType to OrthogPoly types
  template <typename T, typename S>
  struct IsADType< ETPCE::OrthogPoly<T,S> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to OrthogPoly types
  template <typename T, typename S>
  struct IsScalarType< ETPCE::OrthogPoly<T,S> > {
    static const bool value = false;
  };

  //! Specialization of %Value to OrthogPoly types
  template <typename T, typename S>
  struct Value< ETPCE::OrthogPoly<T,S> > {
    typedef typename ValueType< ETPCE::OrthogPoly<T,S> >::type value_type;
    static const value_type& eval(const ETPCE::OrthogPoly<T,S>& x) { 
      return x.val(); }
  };

  //! Specialization of %ScalarValue to OrthogPoly types
  template <typename T, typename S>
  struct ScalarValue< ETPCE::OrthogPoly<T,S> > {
    typedef typename ValueType< ETPCE::OrthogPoly<T,S> >::type value_type;
    typedef typename ScalarType< ETPCE::OrthogPoly<T,S> >::type scalar_type;
    static const scalar_type& eval(const ETPCE::OrthogPoly<T,S>& x) { 
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to OrthogPoly types
  template <typename T, typename S>
  struct StringName< ETPCE::OrthogPoly<T,S> > {
    static std::string eval() { 
      return std::string("Sacado::ETPCE::OrthogPoly< ") + 
	StringName<T>::eval() + " >"; }
  };

  //! Specialization of %IsEqual to OrthogPoly types
  template <typename T, typename S>
  struct IsEqual< ETPCE::OrthogPoly<T,S> > {
    static bool eval(const ETPCE::OrthogPoly<T,S>& x, 
		     const ETPCE::OrthogPoly<T,S>& y) {
      return x.isEqualTo(y);
    }
  };

  //! Specialization of %IsStaticallySized to OrthogPoly types
  template <typename T, typename S>
  struct IsStaticallySized< ETPCE::OrthogPoly<T,S> > {
    static const bool value = false;
  };

} // namespace Sacado

// Define Teuchos traits classes
#ifdef HAVE_SACADO_TEUCHOS
#include "Teuchos_PromotionTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Sacado_PCE_ScalarTraitsImp.hpp"
#include "Teuchos_SerializationTraits.hpp"

namespace Teuchos {

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename T, typename S>
  struct PromotionTraits< Sacado::ETPCE::OrthogPoly<T,S>, 
			  Sacado::ETPCE::OrthogPoly<T,S> > {
    typedef typename Sacado::Promote< Sacado::ETPCE::OrthogPoly<T,S>,
				      Sacado::ETPCE::OrthogPoly<T,S> >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename T, typename S, typename R>
  struct PromotionTraits< Sacado::ETPCE::OrthogPoly<T,S>, R > {
    typedef typename Sacado::Promote< Sacado::ETPCE::OrthogPoly<T,S>, R >::type 
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename L, typename T, typename S>
  struct PromotionTraits< L, Sacado::ETPCE::OrthogPoly<T,S> > {
  public:
    typedef typename Sacado::Promote< L, Sacado::ETPCE::OrthogPoly<T,S> >::type 
    promote;
  };

  //! Specializtion of %Teuchos::ScalarTraits
  template <typename T, typename S>
  struct ScalarTraits< Sacado::ETPCE::OrthogPoly<T,S> > :
    public Sacado::PCE::ScalarTraitsImp< Sacado::ETPCE::OrthogPoly<T,S> > {};

  //! Specialization of %Teuchos::SerializationTraits
  template <typename Ordinal, typename T, typename S>
  struct SerializationTraits<Ordinal, Sacado::ETPCE::OrthogPoly<T,S> > :
    public Sacado::PCE::SerializationTraitsImp< Ordinal, 
						Sacado::ETPCE::OrthogPoly<T,S> >
  {};

  //! Specialization of %Teuchos::ValueTypeSerializer
  template <typename Ordinal, typename T, typename S>
  struct ValueTypeSerializer<Ordinal, Sacado::ETPCE::OrthogPoly<T,S> > :
    public Sacado::PCE::SerializerImp< Ordinal, 
				       Sacado::ETPCE::OrthogPoly<T,S>,
				       ValueTypeSerializer<Ordinal,T> > 
  {
    typedef Sacado::ETPCE::OrthogPoly<T,S> PCEType;
    typedef ValueTypeSerializer<Ordinal,T> ValueSerializer;
    typedef Sacado::PCE::SerializerImp< Ordinal,PCEType,ValueSerializer> Base;
    typedef typename Base::expansion_type expansion_type;
    ValueTypeSerializer(const Teuchos::RCP<expansion_type>& expansion,
			const Teuchos::RCP<const ValueSerializer>& vs) :
      Base(expansion,vs) {}
  };
}
#endif // HAVE_SACADO_TEUCHOS

#endif // SACADO_ETPCE_ORTHOGPOLYTRAITS_HPP
