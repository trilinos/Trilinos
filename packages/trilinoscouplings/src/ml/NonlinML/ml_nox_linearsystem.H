// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
// ML-headers
#include "ml_common.h"
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

#ifndef ML_NOX_LINEARSYSTEM_H
#define ML_NOX_LINEARSYSTEM_H

// NOX
#include "NOX_Epetra_Vector.H"	
#include "NOX_Utils.H"          
#include "NOX_Common.H"         
#include "NOX_EpetraNew_LinearSystem.H"
#include "NOX_EpetraNew_MatrixFree.H"
#include "NOX_EpetraNew_Interface_Jacobian.H"
#include "NOX_EpetraNew_Interface_Preconditioner.H"
#include "NOX_EpetraNew_Interface_Required.H"

// Epetra
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

#include "ml_nox_coarseinterface.H"

// Forward declares
namespace NOX {
  namespace EpetraNew {
    class Scaling;
  }
  namespace Parameter {
    class List;
  }
}


class Epetra_Vector;
class Epetra_Operator;
class Epetra_MultiVector;
class Epetra_Vector;

namespace ML_NOX {

// a linear system class to work with NOX' nlnCG, not a full impl. of a linear system!
class Ml_Nox_LinearSystem : public virtual NOX::EpetraNew::LinearSystem 
{

public:
  //! Constructor.
  Ml_Nox_LinearSystem(NOX::EpetraNew::Interface::Jacobian& iJac,
                      Epetra_Operator& J,
                      NOX::EpetraNew::Interface::Preconditioner& iPrec,
                      ML_NOX::Nox_CoarseProblem_Interface* coarseinterface,
                      Epetra_Operator& Prec,
                      const Epetra_Vector& soln,
                      bool ismatrixfree,
                      int level,
                      int printlevel);

  //! Destructor.
  virtual ~Ml_Nox_LinearSystem();

  /*! 
    \brief Applies Jacobian to the given input vector and puts the answer in the result.

    Computes 
    \f[ v = J u, \f]
    where \f$J\f$ is the Jacobian, \f$u\f$ is the input vector, 
    and \f$v\f$ is the result vector.  Returns true if successful.
  */
  virtual bool applyJacobian(const NOX::Epetra::Vector& input, 
			     NOX::Epetra::Vector& result) const ;

  /*!  
    \brief Applies Jacobian-Transpose to the given input vector and puts the answer in the result.

    Computes 
    \f[ v = J^T u, \f]
    where \f$J\f$ is the Jacobian, \f$u\f$ is the input vector, and \f$v\f$ is the result vector.  Returns true if successful.
    
  */
  virtual bool applyJacobianTranspose(const NOX::Epetra::Vector& input, 
				      NOX::Epetra::Vector& result) const;
  
  /*!
    \brief Applies the inverse of the Jacobian matrix to the given
    input vector and puts the answer in result.
    
    Computes 
    \f[ v = J^{-1} u, \f]
    where \f$J\f$ is the Jacobian, \f$u\f$ is the input vector, 
    and \f$v\f$ is the result vector.
    
    The parameter list contains the linear solver options.
  */
  virtual bool applyJacobianInverse(NOX::Parameter::List &params, 
				    const NOX::Epetra::Vector &input, 
				    NOX::Epetra::Vector &result);

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
				      NOX::Parameter::List& params, 
				      const NOX::Epetra::Vector& input, 
				      NOX::Epetra::Vector& result) const;

  /*! 
    \brief Sets the diagonal scaling vector(s) used in scaling the linear system.  

    See NOX::EpetraNew::Scaling for details on how to specify scaling 
    of the linear system.  
   */
  virtual void resetScaling(NOX::EpetraNew::Scaling& scalingObject);

  //! Evaluates the Jacobian based on the solution vector x.
  virtual bool computeJacobian(Epetra_Vector& x);
 
  /*! 
    \brief Explicitly constructs a preconditioner based on the solution vector x and the parameter list p.  

    The user has the option of recomputing the graph when a new 
    preconditioner is created. The NOX::EpetraNew::Group controls the
    isValid flag for the preconditioner and will control when to call this.
  */
  virtual bool createPreconditioner(Epetra_Vector& x, 
				    NOX::Parameter::List& p,
				    bool recomputeGraph) const;

  /*!
    \brief Deletes the preconditioner. 

    The NOX::EpetraNew::Group controls the isValid flag for the preconditioner and will control when to call this.
  */
  virtual bool destroyPreconditioner() const ;

  /*!
    \brief Allows reuse of a preconditioner, which is handled separately by the ML_Nox_Preconditioner

    The child object will manage setting of this bool flag.  Default implementation is to return false.
  */
  virtual bool checkPreconditionerReuse()
    { return false; } 

  //! Return Jacobian operator
  virtual const Epetra_Operator& getJacobianOperator() const ;

  //! Return Jacobian operator
  virtual Epetra_Operator& getJacobianOperator();

  //! Return preconditioner operator
  virtual const Epetra_Operator& getGeneratedPrecOperator() const;

  //! Return preconditioner operator
  virtual Epetra_Operator& getGeneratedPrecOperator();

  //! Set Jacobian operator for solve
  virtual void 
  setJacobianOperatorForSolve(const Epetra_Operator& solveJacOp);

  //! Set preconditioner operator for solve
  virtual void 
  setPrecOperatorForSolve(const Epetra_Operator& solvePrecOp);

private:

  // the Jacobian interface, referring to either the interface or
  //                         referring to the NOX::EpetraNew::MatrixFree operator that
  //                         implements an NOX::EpetraNew::Interface::Jacobian
  NOX::EpetraNew::Interface::Jacobian& iJac_;              
  
  // the Jacobian, either an Epetra_RowMatrix from the interface
  // or an NOX::EpetraNew::MatrixFree operator
  Epetra_Operator& J_;
  
  // interface referring to the ML_Nox_Preconditioner that implements
  // a NOX::EpetraNew::Interface::Preconditioner and
  // a Epetra_Operator
  NOX::EpetraNew::Interface::Preconditioner& iPrec_;
  
  // the ML_Nox_Preconditioner
  Epetra_Operator*       Precptr_;
  
  // the current solution vector
  const Epetra_Vector& soln_;
  
  // the level this linearsystem lives on
  int level_;
  
  // indicates wether Jacobian is NOX::EpetraNew::MatrixFree operator
  //                           or Epetra_RowMatrix from the interface
  bool ismatrixfree_;
  
  // indicator for output level (0-10)
  int printlevel_;
 
  // a coarse grid interface
  ML_NOX::Nox_CoarseProblem_Interface* coarseinterface_;
  
};

} // namespace ML_NOX


#endif // ML_NOX_LINEARSYSTEM_H
#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA)
