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

#ifndef BELOS_TYPES_HPP
#define BELOS_TYPES_HPP

/*!	\file BelosTypes.hpp
	\brief Collection of the enumerated lists used in Belos.
*/


namespace Belos {

/*! 	\enum Belos::NormType
	\brief Enumerated list for describing the multivector norm type.
*/
  	enum NormType {   OneNorm,       /*!< Compute the one-norm \f$\sum_{i=1}^{n}(|x_i w_i|)\f$ for each vector. */
		    TwoNorm,       /*!< Compute the two-norm *\f$\sqrt(\sum_{i=1}^{n}((x_i w_i)^2)\f$ for each vector. */
		    InfNorm       /*!< Compute the infinity-norm \f$(\max_{i=1}^{n}\{|x_i w_i|\})\f$ for each vector. */
  };

/*!	\enum Belos::ReturnType
	\brief Any method in the Belos abstract interfaces may fail or not be defined. 
	This information needs to be passed back to the algorithm or user.  This will be used 
	by the algorithm or user to decide what should be done.
*/

	enum ReturnType {		Ok, 		/*!< Computation completed sucessfully */
					Undefined, 	/*!< This operation is not defined */
					Error		/*!< This operator returned an error */
	};

/*! \enum Belos::StatusType 
    When the CheckStatus and GetStatus methods of Belos::StatusTest objects are called a 
    variable of type Belos::StatusType is returned.
*/

	enum StatusType { 	Unchecked = 2,   /*!< Initial state of status */
			  	Unconverged = 1, /*!< Convergence is not reached. */
				Converged = 0,   /*!< Convergence is reached. */
				Failed = -1,      /*!< Some failure occured.  Should stop */
			  	NaN = -2         /*!< Result from test contains a NaN value.  Should stop */
			  
	};

} // end Belos namespace

#endif /* BELOS_TYPES_HPP */
