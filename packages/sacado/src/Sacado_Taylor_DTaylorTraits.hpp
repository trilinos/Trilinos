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

#ifndef SACADO_TAYLOR_DTAYLORTRAITS_HPP
#define SACADO_TAYLOR_DTAYLORTRAITS_HPP

#include "Sacado_ConfigDefs.h"
#include "Sacado_ADTraits.hpp"
#include "Sacado_Promote.hpp"

// Forward declarations
namespace Sacado {
  namespace Taylor {
    template <typename T> class DTaylor;
  }
}

namespace Sacado {

  //! Specialization of ADTraits to DTaylor types
  template <typename T>
  class ADTraits< Taylor::DTaylor<T>, Taylor::DTaylor<T> > {
  public:

    typedef Taylor::DTaylor<T> promote;
  };

  //! Specialization of ADTraits to DTaylor types
  template <typename L, typename R>
  class ADTraits< Taylor::DTaylor<L>, R > {
  public:

    typedef typename Taylor::DTaylor<L>::value_type value_type_l;
    typedef typename Promote<R>::value_type value_type_r;
    typedef typename ADTraits<value_type_l,value_type_r>::promote value_type;

    typedef Taylor::DTaylor<value_type> promote;
  };

  //! Specialization of ADTraits to DTaylor types
  template <typename L, typename R>
  class ADTraits< L, Taylor::DTaylor<R> > {
  public:

    typedef typename Promote<L>::value_type value_type_l;
    typedef typename Taylor::DTaylor<R>::value_type value_type_r;
    typedef typename ADTraits<value_type_l,value_type_r>::promote value_type;

    typedef Taylor::DTaylor<value_type> promote;
  };

} // namespace Sacado

#endif // SACADO_TAYLOR_DTAYLORTRAITS_HPP
