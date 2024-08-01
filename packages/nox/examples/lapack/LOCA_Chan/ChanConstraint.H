// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef CHANCONSTRAINT_H
#define CHANCONSTRAINT_H

#include "LOCA.H"
#include "LOCA_LAPACK.H"

class ChanConstraint :
  public LOCA::MultiContinuation::ConstraintInterface {

public:

  //! Constructor
  ChanConstraint(int N, const LOCA::ParameterVector& pVec);

  //! Copy constructor
  ChanConstraint(const ChanConstraint& source,
         NOX::CopyType type = NOX::DeepCopy);

  //! Destructor
  virtual ~ChanConstraint();

  /*!
   * @name Implementation of LOCA::MultiContinuation::ConstraintInterface
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
        NOX::Abstract::MultiVector::DenseMatrix& dgdp,
        bool isValidG);

  //! Return \c true if constraint residuals are valid
  virtual bool isConstraints() const;

  //! Return \c true if derivative of constraint w.r.t. x is valid
  virtual bool isDX() const;

  //! Return constraint residuals
  virtual const NOX::Abstract::MultiVector::DenseMatrix&
  getConstraints() const;

  //! Compute result_p = alpha * dg/dx * input_x
  virtual NOX::Abstract::Group::ReturnType
  multiplyDX(double alpha,
         const NOX::Abstract::MultiVector& input_x,
         NOX::Abstract::MultiVector::DenseMatrix& result_p) const;

  //! Compute result_x = alpha * dg/dx^T * op(b) + beta * result_x
  virtual NOX::Abstract::Group::ReturnType
  addDX(Teuchos::ETransp transb,
    double alpha,
    const NOX::Abstract::MultiVector::DenseMatrix& b,
    double beta,
    NOX::Abstract::MultiVector& result_x) const;

  /*!
   * \brief Return \c true if solution component of constraint
   * derivatives is zero
   */
  virtual bool isDXZero() const;

  //@}

private:

  //! Prohibit generation and use of operator=()
  ChanConstraint& operator=(const ChanConstraint& source);

protected:

  //! Dimension of discretization
  int n;

  //! Constraint values
  NOX::Abstract::MultiVector::DenseMatrix constraints;

  //! Flag indicating whether constraints are valid
  bool isValidConstraints;

  //! Parameter vector
  LOCA::ParameterVector p;

  //! Solution vector
  Teuchos::RCP<NOX::Abstract::MultiVector> x;

};

#endif
