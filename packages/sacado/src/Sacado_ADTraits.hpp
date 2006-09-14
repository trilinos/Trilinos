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

#ifndef SACADO_ADTRAITS_HPP
#define SACADO_ADTRAITS_HPP

#include "Sacado_ConfigDefs.h"

namespace Sacado {

  //! Base template specification for %ADTraits
  /*!
   * The %ADTraits classes provide a mechanism for computing the 
   * promoted type of a binary operation.
   */
  template <typename A, typename B> class ADTraits {};

  //! Specialization of ADTraits for a single type
  template <typename A> class ADTraits<A,A> {
    
  public:
    
    //! Type of promoted values
    typedef A promote;

  };

  //! Specialization of ADTraits to builtin types
#define SACADO_ADTRAITS_SPECIALIZATION(type1,type2,type3) \
  template <> class ADTraits< type1, type2 > {		   \
  public:						   \
    typedef type3 promote;				   \
  };							   \
  template <> class ADTraits< type2, type1 > {		   \
  public:						   \
    typedef type3 promote;				   \
  };

  SACADO_ADTRAITS_SPECIALIZATION(double,float,double)
  SACADO_ADTRAITS_SPECIALIZATION(double,long,double)
  SACADO_ADTRAITS_SPECIALIZATION(double,int,double)
  SACADO_ADTRAITS_SPECIALIZATION(float,long,float)
  SACADO_ADTRAITS_SPECIALIZATION(float,int,float)

} // namespace Sacado

#endif // SACADO_ADTRAITS_HPP
