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

#ifndef SACADO_PROMOTE_HPP
#define SACADO_PROMOTE_HPP

#include "Sacado_ConfigDefs.h"

namespace Sacado {

  //! %Promote class to support Traits
  /*!
   * This class provides a mechanism for all types, AD types and others,
   * to have a consistent ADTraits interface by supplying a \c value_type.
   * Specializations to builtin types are provided.
   */
  template <typename T> class Promote {

  public:

    //! Typename of values
    typedef typename T::value_type value_type;

    //! Constructor
    Promote(const T& x) : x_(x) {}

    //! Return value
    value_type val() { return x_.val(); }

  protected:

    //! AD object
    const T& x_;

  };

  //! Specialization of Promote to a built in type
#define SACADO_PROMOTE_SPECILIZATION(type) \
  template <> class Promote < type >  {	    \
  public:				    \
    typedef type value_type;		    \
    Promote(const type& x) : x_(x) {}	    \
    value_type val() { return x_; }	    \
  protected:				    \
    const type& x_;			    \
  };

  SACADO_PROMOTE_SPECILIZATION(double)
  SACADO_PROMOTE_SPECILIZATION(float)
  SACADO_PROMOTE_SPECILIZATION(long)
  SACADO_PROMOTE_SPECILIZATION(int)

} // namespace Sacado

#endif // SACADO_ADPROMOTE_HPP
