// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_FAD_GENERALFADTRAITS_HPP
#define SACADO_FAD_GENERALFADTRAITS_HPP

#include "Sacado_Traits.hpp"

// Forward declarations
namespace Sacado {
  namespace Fad {
    template <typename T, typename S> class GeneralFad;
  }
}

namespace Sacado {

  //! Specialization of %Promote to GeneralFad types
  template <typename T, typename S>
  struct Promote< Fad::GeneralFad<T,S>, Fad::GeneralFad<T,S> > {
    typedef Fad::GeneralFad<T,S> type;
  };

  //! Specialization of %Promote to GeneralFad types
  template <typename L, typename R, typename S>
  struct Promote< Fad::GeneralFad<L,S>, R > {
    typedef typename ValueType< Fad::GeneralFad<L,S> >::type value_type_l;
    typedef typename Promote<R,R>::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef Fad::GeneralFad<value_type,S> type;
  };

  //! Specialization of %Promote to GeneralFad types
  template <typename L, typename R, typename S>
  struct Promote< L, Fad::GeneralFad<R,S> > {
  public:

    typedef typename Promote<L,L>::type value_type_l;
    typedef typename ValueType< Fad::GeneralFad<R,S> >::type value_type_r;
    typedef typename Promote<value_type_l,value_type_r>::type value_type;

    typedef Fad::GeneralFad<value_type,S> type;
  };

  //! Specialization of %ScalarType to GeneralFad types
  template <typename T, typename S>
  struct ScalarType< Fad::GeneralFad<T,S> > {
    typedef T type;
  };

  //! Specialization of %ValueType to GeneralFad types
  template <typename T, typename S>
  struct ValueType< Fad::GeneralFad<T,S> > {
    typedef T type;
  };

  //! Specialization of %ValueType to GeneralFad types
  template <typename T, typename S>
  struct ScalarValueType< Fad::GeneralFad<T,S> > {
    typedef typename ScalarValueType<T>::type type;
  };

  //! Specialization of %IsADType to GeneralFad types
  template <typename T, typename S>
  struct IsADType< Fad::GeneralFad<T,S> > {
    static const bool value = true;
  };

  //! Specialization of %IsADType to GeneralFad types
  template <typename T, typename S>
  struct IsScalarType< Fad::GeneralFad<T,S> > {
    static const bool value = false;
  };

  //! Specialization of %Value to GeneralFad types
  template <typename T, typename S>
  struct Value< Fad::GeneralFad<T,S> > {
    typedef typename ValueType< Fad::GeneralFad<T,S> >::type value_type;
    static const value_type& eval(const Fad::GeneralFad<T,S>& x) { 
      return x.val(); }
  };

} // namespace Sacado

#endif // SACADO_FAD_GENERALFADTRAITS_HPP
