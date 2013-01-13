// $Id$ 
// $Source$ 

//@HEADER
// ************************************************************************
// 
//            LOCA: Library of Continuation Algorithms Package
//                 Copyright (2005) Sandia Corporation
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

#ifndef LINEARCONSTRAINT_H
#define LINEARCONSTRAINT_H

#include "LOCA.H"

//! A simple linear constraint equation
/*!
 * This class implements the LOCA::MultiContinuation::ConstraintInterfaceMVDX
 * constraint interface for a simple linear constraint:
 * \f[
 *      g(x,p) = \frac{dg}{dx}*x 
 * \f]
 * where \f$x\f$ is the solution vector.
 */
class LinearConstraint : 
  public LOCA::MultiContinuation::ConstraintInterfaceMVDX {
  
public:
  
  //! Constructor
  LinearConstraint(int NumConstraints, const LOCA::ParameterVector& pVec,
		   NOX::Abstract::Vector& cloneVec);

  //! Copy constructor
  LinearConstraint(const LinearConstraint& source, 
		   NOX::CopyType type = NOX::DeepCopy);

  //! Destructor
  virtual ~LinearConstraint();

  /*! 
   * @name Implementation of LOCA::MultiContinuation::ConstraintInterfaceMVDX
   * virtual methods 
   */
  //@{

  //! Copy
  virtual void 
  copy(const LOCA::MultiContinuation::ConstraintInterface& source);

  //! Cloning function
  virtual 
  Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
  clone(NOX::CopyType type = NOX::DeepCopy) const;

  //! Return number of constraints
  virtual int numConstraints() const;

  //! Set the solution vector to y.
  virtual void setX(const NOX::Abstract::Vector& y);

  //! Sets parameter indexed by paramID
  virtual void setParam(int paramID, double val);

  //! Sets parameters indexed by paramIDs
  virtual void setParams(const std::vector<int>& paramIDs, 
			 const NOX::Abstract::MultiVector::DenseMatrix& vals);

  //! Compute continuation constraint equations
  virtual NOX::Abstract::Group::ReturnType
  computeConstraints();

  //! Compute derivative of constraints w.r.t. solution vector x
  virtual NOX::Abstract::Group::ReturnType
  computeDX();

  //! Compute derivative of constraints w.r.t. supplied parameters.
  /*!
   * The first column of \c dgdp should be filled with the constraint
   * residuals \f$g\f$ if \c isValidG is \c false.  If \c isValidG is
   * \c true, then the \c dgdp contains \f$g\f$ on input.
   */
  virtual NOX::Abstract::Group::ReturnType
  computeDP(const std::vector<int>& paramIDs, 
	    NOX::Abstract::MultiVector::DenseMatrix& dp, 
	    bool isValidG);

  //! Return \c true if constraint residuals are valid
  virtual bool isConstraints() const;
  
  //! Return \c true if derivative of constraint w.r.t. x is valid
  virtual bool isDX() const;

  //! Return constraint residuals
  virtual const NOX::Abstract::MultiVector::DenseMatrix&
  getConstraints() const;

  //! Return solution component of constraint derivatives
  virtual const NOX::Abstract::MultiVector*
  getDX() const;

  /*! 
   * \brief Return \c true if solution component of constraint 
   * derivatives is zero
   */
  virtual bool isDXZero() const;

  //@}

  //! Set solution component of constraint derivative
  void setDgDx(const NOX::Abstract::MultiVector& A);

  //! Set parameter component of constraint derivative
  void setDgDp(const NOX::Abstract::MultiVector::DenseMatrix& A);

  /*!
   * \brief Set the flag determinging whether isDXZero() 
   * returns true or false.  This is for testing bordered solvers when B
   * is zero.
   */
  void setIsZeroDX(bool flag);

private:

  //! Prohibit generation and use of operator=()
  LinearConstraint& operator=(const LinearConstraint& source);

protected:

  //! Number of constraints
  int m;

  //! Constraint values
  NOX::Abstract::MultiVector::DenseMatrix constraints;

  //! Flag indicating whether constraints are valid
  bool isValidConstraints;

  //! Solution vector
  Teuchos::RCP<NOX::Abstract::MultiVector> x;

  //! Solution component of constraint derivative
  Teuchos::RCP<NOX::Abstract::MultiVector> dgdx;

  //! Flag indicating whether dgdx is zero
  bool isZeroDgDx;

  //! Parameter vector
  LOCA::ParameterVector p;

  //! Parameter vector as a DenseMatrix (view of p)
  NOX::Abstract::MultiVector::DenseMatrix pvec;

  //! Parameter component of constraint derivative
  NOX::Abstract::MultiVector::DenseMatrix dgdp;

};

#endif
