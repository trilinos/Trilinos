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

#ifndef LOCA_MULTICONTINUATION_ARCLENGTHCONSTRAINT_H
#define LOCA_MULTICONTINUATION_ARCLENGTHCONSTRAINT_H

#include "LOCA_MultiContinuation_ConstraintInterface.H" // base class
#include "NOX_Abstract_MultiVector.H"                   // for dense matrix

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace MultiContinuation {
    class ArcLengthGroup;
  }
}

namespace LOCA { 

  namespace MultiContinuation { 

    /*! 
     * \brief Implementation of LOCA::MultiContinuation::ConstraintInterface
     * for arclength continuation.
     */
    /*!
     * This class implements the arclength constraint equation for
     * pseudo-arclength continuation:
     * \f[
     *     g(x,p,x_0,p_0,x^\ast,p^\ast,v,\Delta s)=
     *         (x-x^\ast)^Tv_x + (p-p^\ast) v_p - \Delta s
     * \f]
     * where \f$v_x\f$, \f$v_p\f$ are the solution and  parameter components 
     * of the predictor direction \f$v\f$ respectively.
     */
    class ArcLengthConstraint : 
      public LOCA::MultiContinuation::ConstraintInterface {

    public:

      //! Constructor
      ArcLengthConstraint(
     const Teuchos::RefCountPtr<LOCA::GlobalData>& global_data,
     const Teuchos::RefCountPtr<LOCA::MultiContinuation::ArcLengthGroup>& grp);

      //! Copy constructor
      ArcLengthConstraint(const ArcLengthConstraint& source, 
			  NOX::CopyType type = NOX::DeepCopy);

      //! Destructor
      ~ArcLengthConstraint();

      //! Assignment
      virtual ArcLengthConstraint&
      operator=(const ArcLengthConstraint& source);

      //! Set pointer to arclength group
      virtual void setArcLengthGroup(const Teuchos::RefCountPtr<LOCA::MultiContinuation::ArcLengthGroup>& grp);

      /*! 
       * @name Implementation of LOCA::MultiContinuation::ConstraintInterface 
       * virtual methods 
       */
      //@{

      virtual ConstraintInterface& 
      operator=(const ConstraintInterface& source);

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
      virtual void setParams(
			 const vector<int>& paramIDs, 
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
      virtual const NOX::Abstract::MultiVector::DenseMatrix&
      getConstraintDerivativesP() const;

      //! Compute result_p = alpha * dg/dx * input_x
      virtual NOX::Abstract::Group::ReturnType
      applyConstraintDerivativesX(
		double alpha, 
		const NOX::Abstract::MultiVector& input_x,
		NOX::Abstract::MultiVector::DenseMatrix& result_p) const;

      //! Compute result_x = alpha * dg/dx^T * op(b) + beta * result_x
      virtual NOX::Abstract::Group::ReturnType
      applyConstraintDerivativesX(
		Teuchos::ETransp transb,
		double alpha, 
		const NOX::Abstract::MultiVector::DenseMatrix& b,
		double beta,
		NOX::Abstract::MultiVector& result_x) const;

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
       * This is currently not defined for arc-length constraints, so we throw
       * an error.
       */
      virtual NOX::Abstract::Group::ReturnType
      computeDgDp(const vector<int>& paramIDs, 
		  NOX::Abstract::MultiVector::DenseMatrix& dgdp, 
		  bool isValidG);

      //@}

    protected:

      //! LOCA global data object
      Teuchos::RefCountPtr<LOCA::GlobalData> globalData;

      //! Pointer to arc-length group
      Teuchos::RefCountPtr<
	LOCA::MultiContinuation::ArcLengthGroup> arcLengthGroup;

      //! Constraint values
      NOX::Abstract::MultiVector::DenseMatrix constraints;

      //! Flag indicating whether constraints are valid
      bool isValidConstraints;

    }; // Class ArcLengthConstraint

  } // namespace MultiContinuation

} // namespace LOCA

#endif // LOCA_MULTICONTINUATION_ARCLENGTHCONSTRAINT_H
