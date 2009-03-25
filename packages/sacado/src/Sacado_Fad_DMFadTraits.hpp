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
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
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
    template <typename T1, typename T2> class DMFad;
  }
}

namespace Sacado {

  //! Specialization of %Promote to DMFad types
  template <typename ValueT, typename ScalarT>
  struct Promote< Fad::DMFad<ValueT,ScalarT>, Fad::DMFad<ValueT,ScalarT> > {
    typedef Fad::DMFad<ValueT,ScalarT> type;
  };

  //! Specialization of %Promote to DMFad types
  template <typename ValueT, typename ScalarT, typename R>
  struct Promote< Fad::DMFad<ValueT,ScalarT>, R > {
    typedef typename ValueType< Fad::DMFad<ValueT,ScalarT> >::type value_type_l;
    typedef typename ValueType<R>::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef Fad::DMFad<value_type,ScalarT> type;
  };

  //! Specialization of %Promote to DMFad types
  template <typename L, typename ValueT, typename ScalarT>
  struct Promote< L, Fad::DMFad<ValueT, ScalarT> > {
  public:

    typedef typename ValueType<L>::type value_type_l;
    typedef typename ValueType< Fad::DMFad<ValueT,ScalarT> >::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef Fad::DMFad<value_type,ScalarT> type;
  };

  //! Specialization of %ScalarType to DMFad types
  template <typename ValueT, typename ScalarT>
  struct ScalarType< Fad::DMFad<ValueT,ScalarT> > {
    typedef ScalarT type;
  };

  //! Specialization of %ValueType to DMFad types
  template <typename ValueT, typename ScalarT>
  struct ValueType< Fad::DMFad<ValueT,ScalarT> > {
    typedef ValueT type;
  };

  //! Specialization of %IsADType to DMFad types
  template <typename ValueT, typename ScalarT>
  struct IsADType< Fad::DMFad<ValueT,ScalarT> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to DMFad types
  template <typename ValueT, typename ScalarT>
  struct IsScalarType< Fad::DMFad<ValueT,ScalarT> > {
    static const bool value = false;
  };

  //! Specialization of %Value to DMFad types
  template <typename ValueT, typename ScalarT>
  struct Value< Fad::DMFad<ValueT,ScalarT> > {
    typedef typename ValueType< Fad::DMFad<ValueT,ScalarT> >::type value_type;
    static const value_type& eval(const Fad::DMFad<ValueT,ScalarT>& x) { 
      return x.val(); }
  };

  //! Specialization of %ScalarValue to DMFad types
  template <typename ValueT, typename ScalarT>
  struct ScalarValue< Fad::DMFad<ValueT,ScalarT> > {
    typedef typename ValueType< Fad::DMFad<ValueT,ScalarT> >::type value_type;
    typedef typename ScalarType< Fad::DMFad<ValueT,ScalarT> >::type scalar_type;
    static const scalar_type& eval(const Fad::DMFad<ValueT,ScalarT>& x) { 
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to DMFad types
  template <typename ValueT, typename ScalarT>
  struct StringName< Fad::DMFad<ValueT,ScalarT> > {
    static std::string eval() { 
      return std::string("Sacado::Fad::DMFad< ") + 
	StringName<ValueT>::eval() + ", " + 
	StringName<ScalarT>::eval() + " >"; }
  };

} // namespace Sacado

// Define Teuchos traits classes
#ifdef HAVE_SACADO_TEUCHOS
#include "Teuchos_PromotionTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Sacado_Fad_ScalarTraitsImp.hpp"

namespace Teuchos {

  //! Specialization of %Teuchos::PromotionTraits to DMFad types
  template <typename ValueT, typename ScalarT>
  struct PromotionTraits< Sacado::Fad::DMFad<ValueT,ScalarT>, 
			  Sacado::Fad::DMFad<ValueT,ScalarT> > {
    typedef typename Sacado::Promote< Sacado::Fad::DMFad<ValueT,ScalarT>,
				      Sacado::Fad::DMFad<ValueT,ScalarT> >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DMFad types
  template <typename ValueT, typename ScalarT, typename R>
  struct PromotionTraits< Sacado::Fad::DMFad<ValueT,ScalarT>, R > {
    typedef typename Sacado::Promote< Sacado::Fad::DMFad<ValueT,ScalarT>,
				      R >::type 
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DMFad types
  template <typename L, typename ValueT, typename ScalarT>
  struct PromotionTraits< L, Sacado::Fad::DMFad<ValueT, ScalarT> > {
  public:
    typedef typename Sacado::Promote< L, 
				      Sacado::Fad::DMFad<ValueT,ScalarT> >::type 
    promote;
  };

  //! Specializtion of Teuchos::ScalarTraits
  template <typename ValueT, typename ScalarT>
  struct ScalarTraits< Sacado::Fad::DMFad<ValueT,ScalarT> > :
    public Sacado::Fad::ScalarTraitsImp< Sacado::Fad::DMFad<ValueT,ScalarT> >
  {};
}
#endif // HAVE_SACADO_TEUCHOS

#endif // SACADO_FAD_DMFADTRAITS_HPP
