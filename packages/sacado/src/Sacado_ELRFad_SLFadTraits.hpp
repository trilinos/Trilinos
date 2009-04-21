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

#ifndef SACADO_ELRFAD_SLFADTRAITS_HPP
#define SACADO_ELRFAD_SLFADTRAITS_HPP

#include "Sacado_Traits.hpp"
#include <sstream>

// Forward declarations
namespace Sacado {
  namespace ELRFad {
    template <typename T, int Num> class SLFad;
  }
}

namespace Sacado {

  //! Specialization of %Promote to SLFad types
  template <typename ValueT, int Num>
  struct Promote< ELRFad::SLFad<ValueT,Num>, 
		  ELRFad::SLFad<ValueT,Num> > {
    typedef ELRFad::SLFad<ValueT,Num> type;
  };

  //! Specialization of %Promote to SLFad types
  template <typename ValueT, int Num, typename R>
  struct Promote< ELRFad::SLFad<ValueT,Num>, R > {
    typedef typename ValueType< ELRFad::SLFad<ValueT,Num> >::type value_type_l;
    typedef typename ValueType<R>::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef ELRFad::SLFad<value_type,Num> type;
  };

  //! Specialization of %Promote to SLFad types
  template <typename L, typename ValueT, int Num>
  struct Promote< L, ELRFad::SLFad<ValueT, Num> > {
  public:

    typedef typename ValueType<L>::type value_type_l;
    typedef typename ValueType< ELRFad::SLFad<ValueT,Num> >::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef ELRFad::SLFad<value_type,Num> type;
  };

  //! Specialization of %ScalarType to SLFad types
  template <typename ValueT, int Num>
  struct ScalarType< ELRFad::SLFad<ValueT,Num> > {
    typedef typename ELRFad::SLFad<ValueT,Num>::ScalarT type;
  };

  //! Specialization of %ValueType to SLFad types
  template <typename ValueT, int Num>
  struct ValueType< ELRFad::SLFad<ValueT,Num> > {
    typedef ValueT type;
  };

  //! Specialization of %IsADType to SLFad types
  template <typename ValueT, int Num>
  struct IsADType< ELRFad::SLFad<ValueT,Num> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to SLFad types
  template <typename ValueT, int Num>
  struct IsScalarType< ELRFad::SLFad<ValueT,Num> > {
    static const bool value = false;
  };

  //! Specialization of %Value to SLFad types
  template <typename ValueT, int Num>
  struct Value< ELRFad::SLFad<ValueT,Num> > {
    typedef typename ValueType< ELRFad::SLFad<ValueT,Num> >::type value_type;
    static const value_type& eval(const ELRFad::SLFad<ValueT,Num>& x) { 
      return x.val(); }
  };

  //! Specialization of %ScalarValue to SLFad types
  template <typename ValueT, int Num>
  struct ScalarValue< ELRFad::SLFad<ValueT,Num> > {
    typedef typename ValueType< ELRFad::SLFad<ValueT,Num> >::type value_type;
    typedef typename ScalarType< ELRFad::SLFad<ValueT,Num> >::type scalar_type;
    static const scalar_type& eval(const ELRFad::SLFad<ValueT,Num>& x) { 
      return ScalarValue<value_type>::eval(x.val()); }
  };

  //! Specialization of %StringName to SLFad types
  template <typename ValueT, int Num>
  struct StringName< ELRFad::SLFad<ValueT,Num> > {
    static std::string eval() { 
      std::stringstream ss;
      ss << "Sacado::ELRFad::SLFad< " 
	 << StringName<ValueT>::eval() << ", " << Num << " >";
      return ss.str(); 
    }
  };

} // namespace Sacado

// Define Teuchos traits classes
#ifdef HAVE_SACADO_TEUCHOS
#include "Teuchos_PromotionTraits.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Sacado_Fad_ScalarTraitsImp.hpp"

namespace Teuchos {

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename ValueT, int Num>
  struct PromotionTraits< Sacado::ELRFad::SLFad<ValueT,Num>, 
			  Sacado::ELRFad::SLFad<ValueT,Num> > {
    typedef typename Sacado::Promote< Sacado::ELRFad::SLFad<ValueT,Num>,
				      Sacado::ELRFad::SLFad<ValueT,Num> >::type
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename ValueT, int Num, typename R>
  struct PromotionTraits< Sacado::ELRFad::SLFad<ValueT,Num>, R > {
    typedef typename Sacado::Promote< Sacado::ELRFad::SLFad<ValueT,Num>,
				      R >::type 
    promote;
  };

  //! Specialization of %Teuchos::PromotionTraits to DFad types
  template <typename L, typename ValueT, int Num>
  struct PromotionTraits< L, Sacado::ELRFad::SLFad<ValueT,Num> > {
  public:
    typedef typename Sacado::Promote< L, 
				      Sacado::ELRFad::SLFad<ValueT,Num> >::type 
    promote;
  };

  //! Specializtion of Teuchos::ScalarTraits
  template <typename ValueT, int Num>
  struct ScalarTraits< Sacado::ELRFad::SLFad<ValueT,Num> > :
    public Sacado::Fad::ScalarTraitsImp< Sacado::ELRFad::SLFad<ValueT,Num> >
  {};
}
#endif // HAVE_SACADO_TEUCHOS

#endif // SACADO_ELRFAD_SFADTRAITS_HPP
