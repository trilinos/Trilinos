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

#ifndef BELOS_NORM_TYPE_HPP
#define BELOS_NORM_TYPE_HPP

/*!	\file BelosNormType.hpp
	\brief Enumerated list for Belos::MultiVec or Belos::Vector norm types.
*/

/*!	\enum Belos::NormType
	\brief The norm method in the Belos::MultiVec or Belos::Vector class may 
	not be defined, this information needs to be passed back to the algorithm.
*/

namespace Belos {

  //! \enum Enumerated list for describing the multivector norm type.
  enum NormType {   OneNorm,       /*!< Compute the one-norm \f$\sum_{i=1}^{n}(|x_i w_i|)\f$ for each vector. */
		    TwoNorm,       /*!< Compute the two-norm *\f$\sqrt(\sum_{i=1}^{n}((x_i w_i)^2)\f$ for each vector. */
		    InfNorm       /*!< Compute the infinity-norm \f$(\max_{i=1}^{n}\{|x_i w_i|\})\f$ for each vector. */
  };

}

#endif
// end of file BelosNormType.hpp
