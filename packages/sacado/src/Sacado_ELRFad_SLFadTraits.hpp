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

// Forward declarations
namespace Sacado {
  namespace ELRFad {
    template <typename T1, int Num, typename T2> class SLFad;
  }
}

namespace Sacado {

  //! Specialization of %Promote to SLFad types
  template <typename ValueT, int Num, typename ScalarT>
  struct Promote< ELRFad::SLFad<ValueT,Num,ScalarT>, 
		  ELRFad::SLFad<ValueT,Num,ScalarT> > {
    typedef ELRFad::SLFad<ValueT,Num,ScalarT> type;
  };

  //! Specialization of %Promote to SLFad types
  template <typename ValueT, int Num, typename ScalarT, typename R>
  struct Promote< ELRFad::SLFad<ValueT,Num,ScalarT>, R > {
    typedef typename ValueType< ELRFad::SLFad<ValueT,Num,ScalarT> >::type value_type_l;
    typedef typename ValueType<R>::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef ELRFad::SLFad<value_type,Num,ScalarT> type;
  };

  //! Specialization of %Promote to SLFad types
  template <typename L, typename ValueT, int Num, typename ScalarT>
  struct Promote< L, ELRFad::SLFad<ValueT, Num, ScalarT> > {
  public:

    typedef typename ValueType<L>::type value_type_l;
    typedef typename ValueType< ELRFad::SLFad<ValueT,Num,ScalarT> >::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef ELRFad::SLFad<value_type,Num,ScalarT> type;
  };

  //! Specialization of %ScalarType to SLFad types
  template <typename ValueT, int Num, typename ScalarT>
  struct ScalarType< ELRFad::SLFad<ValueT,Num,ScalarT> > {
    typedef ScalarT type;
  };

  //! Specialization of %ValueType to SLFad types
  template <typename ValueT, int Num, typename ScalarT>
  struct ValueType< ELRFad::SLFad<ValueT,Num,ScalarT> > {
    typedef ValueT type;
  };

   //! Specialization of %ScalarValueType to SLFad types
  template <typename ValueT, int Num, typename ScalarT>
  struct ScalarValueType< ELRFad::SLFad<ValueT,Num,ScalarT> > {
    typedef typename ScalarValueType< ValueT >::type type;
  };

  //! Specialization of %IsADType to SLFad types
  template <typename ValueT, int Num, typename ScalarT>
  struct IsADType< ELRFad::SLFad<ValueT,Num,ScalarT> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to SLFad types
  template <typename ValueT, int Num, typename ScalarT>
  struct IsScalarType< ELRFad::SLFad<ValueT,Num,ScalarT> > {
    static const bool value = false;
  };

  //! Specialization of %Value to SLFad types
  template <typename ValueT, int Num, typename ScalarT>
  struct Value< ELRFad::SLFad<ValueT,Num,ScalarT> > {
    typedef typename ValueType< ELRFad::SLFad<ValueT,Num,ScalarT> >::type value_type;
    static const value_type& eval(const ELRFad::SLFad<ValueT,Num,ScalarT>& x) { 
      return x.val(); }
  };

} // namespace Sacado

#endif // SACADO_ELRFAD_SFADTRAITS_HPP
