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

#ifndef LOCA_EPETRA_TRANSPOSELINEARSYSTEM_EXPLICITTRANSPOSE_H
#define LOCA_EPETRA_TRANSPOSELINEARSYSTEM_EXPLICITTRANSPOSE_H

#include "LOCA_Epetra_TransposeLinearSystem_AbstractStrategy.H" // base class
#include "EpetraExt_Transpose_RowMatrix.h" // class data element

// forward declarations
class Epetra_RowMatrix;
namespace Teuchos {
  class ParameterList;
}
namespace NOX {
  namespace Epetra {
    class LinearSystem;
    class Scaling;
  }
}
namespace LOCA {
  class GlobalData;
}

namespace LOCA {

  namespace Epetra {

    namespace TransposeLinearSystem {

      /*!
       * \brief Method for solving the transpose of a linear system by
       * explicitly forming the transpose of the matrix.
       */
      /* The class is intialized via the \c solverParams parameter list 
       * argument to the constructor.  The parameters this class recognizes 
       * are:
       * <ul>
       * <li> "Transpose Scaling"  -- 
       *    [Teuchos::RefCountPtr<NOX::Epetra::Scaling>]  - 
       *    Scaling object for transpose solve
       * </ul>
       */
      class ExplicitTranspose : public virtual AbstractStrategy {
	
      public:

	//! Constructor
	ExplicitTranspose(
	    const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
	    const Teuchos::RefCountPtr<Teuchos::ParameterList>& solverParams,
	    const Teuchos::RefCountPtr<NOX::Epetra::LinearSystem>& linsys);

	//! Destructor
	virtual ~ExplicitTranspose();
  
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
				      NOX::Epetra::Vector &result);

  

	//! Evaluates the Jacobian-transpose based on the solution vector x.
	/*!
	 * Note:  For flexibility, this method does not compute the original
	 * Jacobian matrix.  It uses whatever is currently stored in the
	 * linear system.
	 */
	virtual bool 
	createJacobianTranspose();
 
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
				      Teuchos::ParameterList& p);

	//! Get Jacobian-transpose operator
	virtual Teuchos::RefCountPtr<Epetra_Operator> 
	getJacobianTransposeOperator();

	//! Get transpose-preconditioner
	virtual Teuchos::RefCountPtr<Epetra_Operator> 
	getTransposePreconditioner();

	//! Set Jacobian-transpose operator
	virtual void setJacobianTransposeOperator(
	       const Teuchos::RefCountPtr<Epetra_Operator>& new_jac_trans);

	//! Set transpose-preconditioner
	virtual void setTransposePreconditioner(
	      const Teuchos::RefCountPtr<Epetra_Operator>& new_prec_trans);

      protected:

	//! Global data object
	Teuchos::RefCountPtr<LOCA::GlobalData> globalData;

	//! Linear system solver
	Teuchos::RefCountPtr<NOX::Epetra::LinearSystem> linsys;

	//! Jacobian tranpose operator
	Teuchos::RefCountPtr<Epetra_Operator> jac_trans;

	//! Transpose preconditioner
	Teuchos::RefCountPtr<Epetra_Operator> prec_trans;

	//! Scaling object for transpose solve
	Teuchos::RefCountPtr<NOX::Epetra::Scaling> scaling_trans;

	//! Row-matrix transposer
	EpetraExt::RowMatrix_Transpose transposer;

      }; // class TransposePreconditioner

    } // namespace TransposeLinearSystem

  } // namespace Epetra

} // namespace LOCA


#endif
