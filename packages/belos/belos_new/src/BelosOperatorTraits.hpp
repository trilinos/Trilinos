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

#ifndef BELOS_OPERATOR_TRAITS_HPP
#define BELOS_OPERATOR_TRAITS_HPP

/*!     \file BelosOperatorTraits.hpp
        \brief Virtual base class which defines the operator interface 
	required by the iterative linear solver.
*/

#include "BelosMultiVec.hpp"
#include "BelosOperator.hpp"
#include "BelosTypes.hpp"
#include "BelosConfigDefs.hpp"

namespace Belos {

template <class TYPE, class MV, class OP>
class OperatorTraits {};

template <class TYPE> 
class OperatorTraits < TYPE, MultiVec<TYPE>, Operator<TYPE> > 
{
public:

  ///
  static ReturnType Apply ( const Operator<TYPE>& Op, const MultiVec<TYPE>& x, 
			    MultiVec<TYPE>& y, ETrans trans=NOTRANS )
  { return Op.Apply( x, y, trans ); }
  
  ///
  static ReturnType ApplyInverse ( const Operator<TYPE>& Op, const MultiVec<TYPE>& x,
				   MultiVec<TYPE>& y, ETrans trans=NOTRANS )
  { return Op.ApplyInverse( x, y, trans ); }

};
  
  
} // end Belos namespace

#endif // BELOS_OPERATOR_TRAITS_HPP

// end of file BelosOperatorTraits.hpp
