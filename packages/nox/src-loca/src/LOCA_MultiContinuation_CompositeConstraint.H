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

#ifndef LOCA_MULTICONTINUATION_COMPOSITECONSTRAINT_H
#define LOCA_MULTICONTINUATION_COMPOSITECONSTRAINT_H

#include "LOCA_MultiContinuation_ConstraintInterface.H" // base class
#include <vector>

// forward declarations
namespace LOCA {
  class GlobalData;
}

namespace LOCA { 

  namespace MultiContinuation { 

    /*! 
     * \brief Implementation of 
     * LOCA::MultiContinuation::ConstraintInterface for composite
     * constraints, i.e., a constraint comprised of multiple, separate
     * constraints.
     */
    class CompositeConstraint : 
      public virtual LOCA::MultiContinuation::ConstraintInterface {

    public:

      //! Constructor
      CompositeConstraint(
	const Teuchos::RCP<LOCA::GlobalData>& global_data,
	const std::vector< Teuchos::RCP<
	LOCA::MultiContinuation::ConstraintInterface> >& constraintObjects);
      
      //! Copy constructor
      CompositeConstraint(const CompositeConstraint& source, 
			  NOX::CopyType type = NOX::DeepCopy);

      //! Destructor
      ~CompositeConstraint();

      /*! 
       * @name Implementation of LOCA::MultiContinuation::ConstraintInterface 
       * virtual methods 
       */
      //@{

      //! Copy
      virtual void copy(const ConstraintInterface& source);

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
      virtual void setParams(
			 const std::vector<int>& paramIDs, 
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

      //! Return \c true if derivatives of constraints w.r.t. x are valid
      virtual bool isDX() const;

      //! Return constraint residuals
      virtual const NOX::Abstract::MultiVector::DenseMatrix&
      getConstraints() const;

      //! Compute result_p = alpha * dg/dx * input_x
      /*!
       * Note that if there are n constraints and input_x has m columns,
       * result_p should be a n by m matrix and is equivalent to
       * \code
       *    input_x.multiply(alpha, dgdx, result_p).
       * \endcode
       */
      virtual NOX::Abstract::Group::ReturnType
      multiplyDX(double alpha, 
	      const NOX::Abstract::MultiVector& input_x,
	      NOX::Abstract::MultiVector::DenseMatrix& result_p) const;

      //! Compute result_x = alpha * dg/dx^T * op(b) + beta * result_x
      /*!
       * Note that this should be equivalent to
       * \code
       *    result_x.update(transb, alpha, dgdx, b, beta);
       * \endcode
       */
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

      //! Perform any preprocessing before a continuation step starts.
      /*!
       * The \c stepStatus argument indicates whether the previous step was
       * successful.  The default implementation is empty.
       */
      virtual void 
      preProcessContinuationStep(
			     LOCA::Abstract::Iterator::StepStatus stepStatus);
	
      //! Perform any postprocessing after a continuation step finishes.
      /*!
       * The \c stepStatus argument indicates whether the step was
       * successful.  The default implementation is empty.
       */
      virtual void 
      postProcessContinuationStep(
			      LOCA::Abstract::Iterator::StepStatus stepStatus);

      //@}

    private:

      //! Prohibit generation and use of operator=()
      CompositeConstraint& operator=(const CompositeConstraint& source);

    protected:

      //! Default constructor.  Used by derived classes
      CompositeConstraint();

      //! Initialize object when using default constructor
     virtual void init(
	const Teuchos::RCP<LOCA::GlobalData>& global_data,
	const std::vector< Teuchos::RCP<
	LOCA::MultiContinuation::ConstraintInterface> >& constraintObjects);

    protected:

      //! LOCA global data object
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Number of constraint objects
      int numConstraintObjects;

      //! Array of constraint pointers
      std::vector< Teuchos::RCP<
	LOCA::MultiContinuation::ConstraintInterface> > constraintPtrs;

      //! Indices of each sub-constraint in the composite constraint
      std::vector< std::vector<int> > indices;

      //! Total number of constraints
      int totalNumConstraints;
      
      //! Constraint values
      NOX::Abstract::MultiVector::DenseMatrix constraints;

      //! Flag indicating whether constraints are valid
      bool isValidConstraints;

      //! Flag indicating whether dg/dx is valid
      bool isValidDX;

    }; // Class CompositeConstraint

  } // namespace MultiContinuation

} // namespace LOCA

#endif // LOCA_MULTICONTINUATION_COMPOSITECONSTRAINT_H
