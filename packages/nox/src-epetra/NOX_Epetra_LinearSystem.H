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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

#ifndef NOX_EPETRA_LINEARSYSTEM_H
#define NOX_EPETRA_LINEARSYSTEM_H

#include "NOX_Epetra_Vector.H"	    // class data element
#include "NOX_Utils.H"              // class data element
#include "NOX_Common.H"             // class data element (std::string)
#include "Teuchos_RCP.hpp"  // class data element

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
  //! Determines handling of the preconditioner between nonlinear iterations.
  enum PreconditionerReusePolicyType {
    //! Destroy and recreate the preconditioner between nonlinear iterations.
    PRPT_REBUILD,
    //! Recompute using already allocated structures for preconditioner.
    PRPT_RECOMPUTE,
    //! Reuse the preconditioner from previous iteration.
    PRPT_REUSE
  };  

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
  virtual Teuchos::RCP<NOX::Epetra::Scaling> getScaling() = 0;

  /*! 
    \brief Sets the diagonal scaling vector(s) used in scaling the linear system.  

    See NOX::Epetra::Scaling for details on how to specify scaling 
    of the linear system.  
   */
  virtual void 
  resetScaling(const Teuchos::RCP<NOX::Epetra::Scaling>& s) = 0;

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

  /*! \brief Recalculates the preconditioner using an already allocated graph.
    
    Use this to compute a new preconditioner while using the same
    graph for the preconditioner.  This avoids deleting and
    reallocating the memory required for the preconditioner and
    results in a big speed-up for large-scale jobs.
  */
  virtual bool recomputePreconditioner(const NOX::Epetra::Vector& x, 
			Teuchos::ParameterList& linearSolverParams) const = 0;

  /*! \brief  Evaluates the preconditioner policy at the current state.

     NOTE: This can change values between nonlienar iterations.  It is
     not a static value.
  */
  virtual PreconditionerReusePolicyType 
  getPreconditionerPolicy(bool advanceReuseCounter=true) = 0;

  //! Indicates whether a preconditioner has been constructed
  virtual bool isPreconditionerConstructed() const = 0;

  //! Indicates whether the linear system has a preconditioner
  virtual bool hasPreconditioner() const = 0;

  //! Return Jacobian operator
  virtual Teuchos::RCP<const Epetra_Operator> 
  getJacobianOperator() const = 0;

  //! Return Jacobian operator
  virtual Teuchos::RCP<Epetra_Operator> getJacobianOperator() = 0;

  //! Return preconditioner operator
  /*!
   * Note:  This should only be called if hasPreconditioner() returns true.
   */
  virtual Teuchos::RCP<const Epetra_Operator> 
  getGeneratedPrecOperator() const = 0;

  //! Return preconditioner operator
  virtual Teuchos::RCP<Epetra_Operator> getGeneratedPrecOperator() = 0;

  //! Set Jacobian operator for solve
  virtual void setJacobianOperatorForSolve(const 
		 Teuchos::RCP<const Epetra_Operator>& solveJacOp) = 0;

  //! Set preconditioner operator for solve
  /*!
   * Note:  This should only be called if hasPreconditioner() returns true.
   */
  virtual void setPrecOperatorForSolve(const 
		 Teuchos::RCP<const Epetra_Operator>& solvePrecOp) = 0;

  //! Statistics for number of times the linear solver has been called (def: 0)
  virtual int getNumLinearSolves() {return 0;};

  //! Statistics for number of iterations taken in last linear solve (def: 0)
  virtual int getLinearItersLastSolve() {return 0;};

  //! Statistics for cumulative number of iterations in all linear solve (def: 0)
  virtual int getLinearItersTotal() {return 0;};

  //! Statistics for the achieved tolerance of last linear solve (def: 0.0)
  virtual double getAchievedTol() {return 0.0;};
};

} // namespace Epetra
} // namespace NOX


#endif
