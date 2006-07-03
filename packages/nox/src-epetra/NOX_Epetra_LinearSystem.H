// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Tammy Kolda (tgkolda@sandia.gov) or Roger Pawlowski
// (rppawlo@sandia.gov).
// 
// ************************************************************************
//@HEADER

#ifndef NOX_EPETRA_LINEARSYSTEM_H
#define NOX_EPETRA_LINEARSYSTEM_H

#include "NOX_Epetra_Vector.H"	    // class data element
#include "NOX_Utils.H"              // class data element
#include "NOX_Common.H"             // class data element (string)
#include "Teuchos_RefCountPtr.hpp"  // class data element

// Forward declares
namespace NOX {
  namespace Epetra {
    class Scaling;
  }
  namespace Parameter {
    class List;
  }
}
class Epetra_Vector;
class Epetra_Operator;

namespace NOX {
namespace Epetra {

//! Pure virtual class interface for allowing different linear solvers to be used by the NOX::Epetra::Group.  
class LinearSystem {

public:
  //! Constructor.
  LinearSystem(){};

  //! Destructor.
  virtual ~LinearSystem(){};

  /*! 
    \brief Applies Jacobian to the given input vector and puts the answer in the result.

    Computes 
    \f[ v = J u, \f]
    where \f$J\f$ is the Jacobian, \f$u\f$ is the input vector, 
    and \f$v\f$ is the result vector.  Returns true if successful.
  */
  virtual bool applyJacobian(const NOX::Epetra::Vector& input, 
			     NOX::Epetra::Vector& result) const = 0;

  /*!  
    \brief Applies Jacobian-Transpose to the given input vector and puts the answer in the result.

    Computes 
    \f[ v = J^T u, \f]
    where \f$J\f$ is the Jacobian, \f$u\f$ is the input vector, and \f$v\f$ is the result vector.  Returns true if successful.
    
  */
  virtual bool applyJacobianTranspose(const NOX::Epetra::Vector& input, 
				      NOX::Epetra::Vector& result) const = 0;
  
  /*!
    \brief Applies the inverse of the Jacobian matrix to the given
    input vector and puts the answer in result.
    
    Computes 
    \f[ v = J^{-1} u, \f]
    where \f$J\f$ is the Jacobian, \f$u\f$ is the input vector, 
    and \f$v\f$ is the result vector.
    
    The parameter list contains the linear solver options.
  */
  virtual bool applyJacobianInverse(Teuchos::ParameterList &params, 
				    const NOX::Epetra::Vector &input, 
				    NOX::Epetra::Vector &result) = 0;

  /*!
    \brief Apply right preconditiong to the given input vector

    Let \f$M\f$ be a right preconditioner for the Jacobian \f$J\f$; in
    other words, \f$M\f$ is a matrix such that
    \f[ JM \approx I. \f]
    
    Compute
    \f[ u = M^{-1} v, \f]
    where \f$u\f$ is the input vector and \f$v\f$ is the result vector.
    
    If <em>useTranspose</em> is true, then the transpose of the 
    preconditioner is applied:  
    \f[ u = {M^{-1}}^T v, \f]
    The transpose preconditioner is currently only required for 
    Tensor methods.

    The parameter list contains the linear solver options.
  */
  virtual bool applyRightPreconditioning(bool useTranspose,
				      Teuchos::ParameterList& params, 
				      const NOX::Epetra::Vector& input, 
				      NOX::Epetra::Vector& result) const = 0;

  //! Get the scaling object
  virtual Teuchos::RefCountPtr<NOX::Epetra::Scaling> getScaling() = 0;

  /*! 
    \brief Sets the diagonal scaling vector(s) used in scaling the linear system.  

    See NOX::Epetra::Scaling for details on how to specify scaling 
    of the linear system.  
   */
  virtual void 
  resetScaling(const Teuchos::RefCountPtr<NOX::Epetra::Scaling>& s) = 0;

  //! Evaluates the Jacobian based on the solution vector x.
  virtual bool computeJacobian(const NOX::Epetra::Vector& x) = 0;
 
  /*! 
    \brief Explicitly constructs a preconditioner based on the solution vector x and the parameter list p.  

    The user has the option of recomputing the graph when a new 
    preconditioner is created. The NOX::Epetra::Group controls the
    isValid flag for the preconditioner and will control when to call this.
  */
  virtual bool createPreconditioner(const NOX::Epetra::Vector& x, 
				    Teuchos::ParameterList& p,
				    bool recomputeGraph) const = 0;

  /*!
    \brief Deletes the preconditioner. 

    The NOX::Epetra::Group controls the isValid flag for the preconditioner and will control when to call this.
  */
  virtual bool destroyPreconditioner() const = 0;

  /*!
    \brief Allows reuse of a preconditioner

    The child object will manage setting of this bool flag.  Default implementation is to return false.
  */
  virtual bool checkPreconditionerReuse()
    { return false; } 

  //! Indicates whether a preconditioner has been constructed
  virtual bool isPreconditionerConstructed() const = 0;

  //! Indicates whether the linear system has a preconditioner
  virtual bool hasPreconditioner() const = 0;

  //! Return Jacobian operator
  virtual Teuchos::RefCountPtr<const Epetra_Operator> 
  getJacobianOperator() const = 0;

  //! Return Jacobian operator
  virtual Teuchos::RefCountPtr<Epetra_Operator> getJacobianOperator() = 0;

  //! Return preconditioner operator
  /*!
   * Note:  This should only be called if hasPreconditioner() returns true.
   */
  virtual Teuchos::RefCountPtr<const Epetra_Operator> 
  getGeneratedPrecOperator() const = 0;

  //! Return preconditioner operator
  virtual Teuchos::RefCountPtr<Epetra_Operator> getGeneratedPrecOperator() = 0;

  //! Set Jacobian operator for solve
  virtual void setJacobianOperatorForSolve(const 
		 Teuchos::RefCountPtr<const Epetra_Operator>& solveJacOp) = 0;

  //! Set preconditioner operator for solve
  /*!
   * Note:  This should only be called if hasPreconditioner() returns true.
   */
  virtual void setPrecOperatorForSolve(const 
		 Teuchos::RefCountPtr<const Epetra_Operator>& solvePrecOp) = 0;
};

} // namespace Epetra
} // namespace NOX


#endif
