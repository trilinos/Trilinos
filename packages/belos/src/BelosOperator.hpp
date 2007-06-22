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

#ifndef BELOS_OPERATOR_HPP
#define BELOS_OPERATOR_HPP

/*!     \file BelosOperator.hpp
        \brief Virtual base class which defines the operator interface 
	required by the iterative linear solver.
*/

#include "BelosOperatorTraits.hpp"
#include "BelosMultiVec.hpp"
#include "BelosConfigDefs.hpp"

/*!	\class Belos::Operator

	\brief Belos's templated pure virtual class for constructing the operator that is
	used by the linear solver.  

	This operator is used as the interface to the matrix (<tt>A</tt>), 
	solution (<tt>X</tt>), and right-hand side (<tt>B</tt>) of the linear system <tt>AX = B</tt>.
	Furthermore, it is also the interface to left/right preconditioning and left/right scaling of the
	linear system.

	A concrete implementation of this class is necessary.  The user can create their own implementation
	if those supplied are not suitable for their needs.

	\author Michael Heroux and Heidi Thornquist
*/

namespace Belos {
  
  template <class ScalarType>
  class Operator {
  public:
    
    //! @name Constructor/Destructor
    //@{ 
    
    //! Default constructor
    Operator() {};
    
    //! Destructor.
    virtual ~Operator() {};
    //@}
    
    //! @name Operator application method
    //@{ 
    
    /*! \brief This routine takes the Belos::MultiVec \c x and applies the operator
      to it resulting in the Belos::MultiVec \c y, which is returned.
        \note It is expected that any problem with applying this operator to \c x will be
	indicated by an exception being thrown.
    */
    virtual void Apply ( const MultiVec<ScalarType>& x, 
			 MultiVec<ScalarType>& y, ETrans trans=NOTRANS ) const = 0;
  };
  
  ////////////////////////////////////////////////////////////////////
  //
  // Implementation of the Belos::OperatorTraits for Belos::Operator 
  //                                               and Belos::MultiVec.
  //
  ////////////////////////////////////////////////////////////////////  
  
  /*!  \brief Template specialization of Belos::OperatorTraits class using Belos::Operator and Belos::MultiVec virtual
    base classes.
    
    Any class that inherits from Belos::Operator will be accepted by the Belos templated solvers due to this
    interface to the Belos::OperatorTraits class.
  */

  template <class ScalarType> 
  class OperatorTraits < ScalarType, MultiVec<ScalarType>, Operator<ScalarType> > 
  {
  public:
    
    ///
    static void Apply ( const Operator<ScalarType>& Op, 
			const MultiVec<ScalarType>& x, 
			MultiVec<ScalarType>& y,
			ETrans trans=NOTRANS )
    { Op.Apply( x, y, trans ); }
    
  };
  
} // end Belos namespace

#endif

// end of file BelosOperator.hpp
