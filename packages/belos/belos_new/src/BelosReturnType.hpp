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

#ifndef BELOS_RETURN_TYPE_HPP
#define BELOS_RETURN_TYPE_HPP

/*!	\file BelosReturnType.hpp
	\brief Enumerated list for Belos::Operator return types.
*/

/*!	\enum Belos::ReturnType
	\brief The apply method in the Belos::Operator may fail or not be defined, 
	this information needs to be passed back to the algorithm.  This will be used 
	in the algorithm to decide what should be done if the Belos::Operator fails for
	any reason.
*/

namespace Belos {

	enum ReturnType {		Ok, 		/*!< Computation completed sucessfully */
					Undefined, 	/*!< This operation is not defined */
					Error		/*!< This operator returned an error */
	};

}

#endif
// end of file BelosReturnType.hpp
