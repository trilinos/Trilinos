// @HEADER
// ***********************************************************************
//
//                 Belos: Block Linear Solvers Package
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

#ifndef ANASAZI_PRECONDITION_HPP
#define ANASAZI_PRECONDITION_HPP

#include "AnasaziMultiVec.hpp"
#include "BelosConfigDefs.hpp"

/*!	\class Anasazi::Precondition

	\brief Belos's templated virtual class for constructing the 
	preconditioner to the AnasaziMatrix that is used by the linear 
	solver.

	A concrete implementation of this class is not necessary.  The user 
	can create their own implementation if those supplied are not
	suitable for their needs.  By default, this class will provide
	the identity as the preconditioner.  Thus, if the ApplyPrecondition
	function is not overridden, then the block linear solver will be
	unpreconditioned.

	\author Teri Barth, Mike Heroux, and Heidi Thornquist
*/

namespace Anasazi {

template <class TYPE>
class Precondition {
public:
	//@{ \name Constructor/Destructor.

	//! %Anasazi::Precondition constructor.
	Precondition() {};

	//! %Anasazi::Precondition destructor.
	virtual ~Precondition() {};
	//@}

	//@{ \name Preconditioner application method.
	
	/*! \brief This routine takes the %Anasazi::MultiVec \c x and 
	applies the preconditioner to it resulting in the %Anasazi::MultiVec \c y, 
	which is returned.  If this routine is not overridden, then the %Anasazi::MultiVec
	\c x will be passed directly to \c y, resulting in an unpreconditioned linear
	solver.
	*/
	virtual void ApplyPrecondition (const MultiVec<TYPE>& x, MultiVec<TYPE>& y ) const {
		//
		// First cast away the const on x.
		//
		MultiVec<TYPE>& temp_x = const_cast<MultiVec<TYPE>& >(x);
		//
		// Now create the indexing for copying x into y.
		//
		int i, numvecs = x.GetNumberVecs();
		int *index = new int[numvecs]; assert(index!=NULL);
		for (i=0; i<numvecs; i++) { index[i] = i; }
		y.SetBlock(temp_x, index, numvecs);
		delete [] index;
	}
	//@}
};

}
#endif
// end of file AnasaziPrecondition.hpp
