// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
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

#ifndef ANASAZI_OPERATOR_HPP
#define ANASAZI_OPERATOR_HPP

#include "AnasaziMultiVec.hpp"
#include "AnasaziReturnType.hpp"
#include "AnasaziConfigDefs.hpp"
#include "Teuchos_ScalarTraits.hpp"

/*!	\class Anasazi::Operator

	\brief Anasazi's templated virtual class for constructing the operator that is
	used by the eigensolver.
	
	A concrete implementation of this class is necessary.  The user can create their own implementation
	if those supplied are not suitable for their needs.

	\author Rich Lehoucq, Heidi Thornquist
*/

namespace Anasazi {

template <class TYPE>
class Operator {
public:
	//@{ \name Constructor/Destructor.
	//! Default constructor.
	Operator() {};

	//! Destructor.
	virtual ~Operator(void) {};
	//@}

	//@{ \name Matrix/Operator application method.

	/*! \brief This routine takes the %Anasazi::MultiVec \c x and
	applies the operator to it resulting in the %Anasazi::MultiVec \c y,
	which is returned.  If this routine is not overridden, then the %Anasazi::MultiVec
	\c x will be passed directly to \c y.  Thus the operator is the identity if this
	method is defined by the user.
	*/
	virtual ReturnType Apply (const MultiVec<TYPE>& x, MultiVec<TYPE>& y ) const 
	{
            if (x.GetNumberVecs() == y.GetNumberVecs()) {	
                //
                // First cast away the const on x.
                //
                MultiVec<TYPE>& temp_x = const_cast<MultiVec<TYPE>& >(x);
                //
                // Now create the indexing for copying x into y.
                //
		TYPE one = Teuchos::ScalarTraits<TYPE>::one();
		TYPE zero = Teuchos::ScalarTraits<TYPE>::zero();
		y.MvAddMv( one, temp_x, zero, temp_x );
		return Ok;
	    }
	    else { return Failed; }
        };
	//@}
};

} // end of Anasazi namespace
#endif
// end of file AnasaziOperator.hpp
