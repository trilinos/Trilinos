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

#ifndef LOCA_EPETRA_TRANSPOSELINEARSYSTEM_ABSTRACTSTRATEGY_H
#define LOCA_EPETRA_TRANSPOSELINEARSYSTEM_ABSTRACTSTRATEGY_H

#include "Teuchos_RefCountPtr.hpp"

// Forward declares
class Epetra_Operator;
namespace Teuchos {
  class ParameterList;
}
namespace NOX {
  namespace Epetra {
    class Vector;
  }
}

namespace LOCA {

  namespace Epetra {

    //! Strategies for solving transpose linear systems
    namespace TransposeLinearSystem {

      //! A pure virtual interface for solving the transpose of a linear system
      class AbstractStrategy {
	
      public:

	//! Constructor
	AbstractStrategy() {}

	//! Destructor
	virtual ~AbstractStrategy() {}
  
	/*!
	 * \brief Applies the inverse of the Jacobian matrix transpose to the 
	 * given input vector and puts the answer in result.
	 *
	 * Computes 
	 * \f[ v = J^{-T} u, \f]
	 * where \f$J\f$ is the Jacobian, \f$u\f$ is the input vector, 
	 * and \f$v\f$ is the result vector.
	 *
	 * The parameter list contains the linear solver options.
	 */
	virtual bool 
	applyJacobianTransposeInverse(Teuchos::ParameterList &params, 
				      const NOX::Epetra::Vector &input, 
				      NOX::Epetra::Vector &result) = 0;

  

	//! Evaluates the Jacobian-transpose based on the solution vector x.
	/*!
	 * Note:  For flexibility, this method does not compute the original
	 * Jacobian matrix.  It uses whatever is currently stored in the
	 * linear system.
	 */
	virtual bool 
	createJacobianTranspose() = 0;
 
	/*! 
	 * \brief Explicitly constructs a preconditioner based on the solution 
	 * vector x and the parameter list p.  
	 */
	/*!
	 * Note:  \c x is only needed for user-supplied preconditioners.
	 * When using a built-in preconditioner (e.g., Ifpack), \c x will
	 * note be used.
	 */
	virtual bool 
	createTransposePreconditioner(const NOX::Epetra::Vector& x, 
				      Teuchos::ParameterList& p) = 0;

	//! Get Jacobian-transpose operator
	virtual Teuchos::RefCountPtr<Epetra_Operator> 
	getJacobianTransposeOperator() = 0;

	//! Get transpose-preconditioner
	virtual Teuchos::RefCountPtr<Epetra_Operator> 
	getTransposePreconditioner() = 0;

	//! Set Jacobian-transpose operator
	virtual void setJacobianTransposeOperator(
	       const Teuchos::RefCountPtr<Epetra_Operator>& new_jac_trans) = 0;

	//! Set transpose-preconditioner
	virtual void setTransposePreconditioner(
	      const Teuchos::RefCountPtr<Epetra_Operator>& new_prec_trans) = 0;

      }; // class AbstractStrategy

    } // namespace TransposeLinearSystem

  } // namespace Epetra

} // namespace LOCA


#endif
