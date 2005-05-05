// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
//
//                  LOCA Continuation Algorithm Package
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
// Questions? Contact Andy Salinger (agsalin@sandia.gov) or Eric Phipps
// (etphipp@sandia.gov), Sandia National Laboratories.
//
// ************************************************************************
//@HEADER

#ifndef NORMCONSTRAINT_H
#define NORMCONSTRAINT_H

#include "LOCA.H"
#include "LOCA_LAPACK.H"

//! A simple constraint equation based on the solution norm
/*!
 * This class implements the LOCA::MultiContinuation::ConstraintInterfaceMVDX
 * constraint interface for a simple norm constraint:
 * \f[
 *      g(x,p) = \frac{1}{2 N} \|x\|^2 - \gamma
 * \f]
 * where \f$x\f$ is the solution vector and \gamma is the constraint
 * parameter.
 */
class NormConstraint : 
  public LOCA::MultiContinuation::ConstraintInterfaceMVDX {
  
public:
  
  //! Constructor
  NormConstraint(int N, const LOCA::ParameterVector& pVec);

  //! Copy constructor
  NormConstraint(const NormConstraint& source, 
		 NOX::CopyType type = NOX::DeepCopy);

  //! Destructor
  virtual ~NormConstraint();

  //! Assignment
  virtual NormConstraint& operator=(const NormConstraint& source);

  /*! 
   * @name Implementation of LOCA::MultiContinuation::ConstraintInterfaceMVDX
   * virtual methods 
   */
  //@{

  virtual LOCA::MultiContinuation::ConstraintInterface& 
  operator=(const LOCA::MultiContinuation::ConstraintInterface& source);

  //! Cloning function
  virtual 
  Teuchos::RefCountPtr<LOCA::MultiContinuation::ConstraintInterface>
  clone(NOX::CopyType type = NOX::DeepCopy) const;

  //! Return number of constraints
  virtual int numConstraints() const;

  //! Set the solution vector to y.
  virtual void setX(const NOX::Abstract::Vector& y);

  //! Sets parameter indexed by paramID
  virtual void setParam(int paramID, double val);

  //! Sets parameters indexed by paramIDs
  virtual void setParams(const vector<int>& paramIDs, 
			 const NOX::Abstract::MultiVector::DenseMatrix& vals);

  //! Compute continuation constraint equations
  virtual NOX::Abstract::Group::ReturnType
  computeConstraints();

  //! Compute continuation constraint equations
  virtual NOX::Abstract::Group::ReturnType
  computeConstraintDerivatives();

  //! Return \c true if constraint residuals are valid
  virtual bool isConstraints() const;
  
  //! Return \c true if constraint derivatives are valid
  virtual bool isConstraintDerivatives() const;

  //! Return constraint residuals
  virtual const NOX::Abstract::MultiVector::DenseMatrix&
  getConstraints() const;

  //! Return parameter component of constraint derivatives
  virtual const NOX::Abstract::MultiVector::DenseMatrix*
  getConstraintDerivativesP() const;

  //! Return solution component of constraint derivatives
  virtual const NOX::Abstract::MultiVector*
  getConstraintDerivativesX() const;

  /*! 
   * \brief Return \c true if solution component of constraint 
   * derivatives is zero
   */
  virtual bool isConstraintDerivativesXZero() const;
  
  /*! 
   * \brief Return \c true if parameter component of constraint 
   * derivatives is zero
   */
  virtual bool isConstraintDerivativesPZero() const;

  //! Compute derivative of constraints w.r.t. supplied parameters.
  /*!
   * This is currently not defined for natural constraints, so we throw
   * an error.
   */
  virtual NOX::Abstract::Group::ReturnType
  computeDgDp(const vector<int>& paramIDs, 
	      NOX::Abstract::MultiVector::DenseMatrix& dgdp, 
	      bool isValidG);

  //@}

protected:

  //! Dimension of discretization
  int n;

  //! Constraint values
  NOX::Abstract::MultiVector::DenseMatrix constraints;

  //! Flag indicating whether constraints are valid
  bool isValidConstraints;

  //! Derivative of constraints w.r.t. constraint parameters
  NOX::Abstract::MultiVector::DenseMatrix dgdp;

  //! Parameter vector
  LOCA::ParameterVector p;

  //! Solution vector
  Teuchos::RefCountPtr<NOX::Abstract::MultiVector> x;

};

#endif
