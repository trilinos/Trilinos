// $Id$
// $Source$

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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
// 
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or 
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

#ifndef LOCA_BORDEREDSOLVER_ABSTRACTOPERATOR_H
#define LOCA_BORDEREDSOLVER_ABSTRACTOPERATOR_H

#include "NOX_Abstract_Group.H"           // for ReturnType
#include "NOX_Abstract_MultiVector.H"

namespace LOCA {

  namespace BorderedSolver {
    
    /*! 
     * \brief Abstract interface class representing an operator for solving 
     * bordered sets of linear equations */
    class AbstractOperator {

    public:

      //! Constructor
      AbstractOperator() {}

      //! Destructor
      virtual ~AbstractOperator() {}
  
      //! Apply the operator
      virtual NOX::Abstract::Group::ReturnType 
      apply(const NOX::Abstract::MultiVector& X,
	    NOX::Abstract::MultiVector& Y) const = 0;

      //! Apply transpose of the operator
      virtual NOX::Abstract::Group::ReturnType 
      applyTranspose(const NOX::Abstract::MultiVector& X,
		     NOX::Abstract::MultiVector& Y) const = 0;

      //! Apply inverse of the operator
      virtual NOX::Abstract::Group::ReturnType 
      applyInverse(Teuchos::ParameterList& params,
		   const NOX::Abstract::MultiVector& B,
		   NOX::Abstract::MultiVector& X) const = 0;

      //! Apply inverse transpose of the operator
      virtual NOX::Abstract::Group::ReturnType 
      applyInverseTranspose(Teuchos::ParameterList& params,
			    const NOX::Abstract::MultiVector& B,
			    NOX::Abstract::MultiVector& X) const = 0;

    private:

      //! Private to prohibit copying
      AbstractOperator(const AbstractOperator&);

      //! Private to prohibit copying
      AbstractOperator& operator = (const AbstractOperator&);

    }; // Class AbstractOperator

  } // Namespace BorderedSolver

} // Namespace LOCA

#endif // LOCA_BORDEREDSOLVER_ABSTRACTOPERATOR_H
