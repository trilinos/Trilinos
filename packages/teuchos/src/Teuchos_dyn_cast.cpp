// @HEADER
// ***********************************************************************
// 
//                    Teuchos: Common Tools Package
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#include "Teuchos_dyn_cast.hpp"
#include "Teuchos_TestForException.hpp"

/** We throw a m_bad_cast, which is a subclass of bad_cast.  
	This is necessary, since bad_cast lacks the appropriate
	constructor for use with the TEST_FOR_EXCEPTION macro.
*/
void Teuchos::dyn_cast_throw_exception(
	const char T_from[], const char T_from_concr[], const char T_to[]
  )
{
	TEST_FOR_EXCEPTION(
		true, m_bad_cast
		,"dyn_cast<" << T_to << ">(" << T_from
		<< ") : Error, the object with the concrete type \'"
		<< T_from_concr << "\' (passed in through the interface type \'" << T_from <<  "\') "
    " does not support the interface \'"
		<< T_to << "\' and the dynamic cast failed!" );
}
