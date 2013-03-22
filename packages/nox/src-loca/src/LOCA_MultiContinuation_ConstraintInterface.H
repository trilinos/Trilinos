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

#ifndef LOCA_MULTICONTINUATION_CONSTRAINTINTERFACE_H
#define LOCA_MULTICONTINUATION_CONSTRAINTINTERFACE_H

#include "Teuchos_RCP.hpp"           // for RCP

#include "NOX_Abstract_Group.H"              // for CopyType, ReturnType
#include "NOX_Abstract_MultiVector.H"        // for MultiVector, DenseMatrix
#include "LOCA_Abstract_Iterator.H"          // for StepStatus

// forward declarations
namespace LOCA {
  namespace MultiContinuation {
    class ExtendedVector;
  }
}

namespace LOCA {

  namespace MultiContinuation {

    /*! 
     * \brief Abstract interface for the constraint portion of a constrained 
     * nonlinear system.
     */
    /*!
     * This class is used in conjunction with 
     * LOCA::MultiContinuation::ConstrainedGroup to represent a constrained
     * nonlinear system:
     * \f[
     *      f(x,y) = 0
     *      g(x,y) = 0
     * \f]
     * where \f$f(x,y)\f$ is represented by a concrete implementation of
     * a LOCA::MultiContinuation::AbstractGroup and \f$g(x,y)\f$ (the
     * constraint) is represented by an implementation of this class.  Here
     * it is assumed the resulting system is square, i.e., \f$x\in\Re^n\f$,
     * \f$y\in\Re^m\f$, \f$f(x,y)\in\Re^n\f$ and \f$g(x,y)\in\Re^m\f$.
     *
     * This class provides an interface to evaluate \f$g(x,y)\f$, compute
     * the derivatives \f$g_x\f$ and \f$g_y\f$, and apply the derivative
     * \f$g_x\f$ to arbitrary multi-vectors (the implementation is never
     * required to explicitly store \f$g_x\f$ which is impractical in 
     * many situations).
     */
    class ConstraintInterface {

    public:

      //! Constructor
      ConstraintInterface() {}

      //! Destructor
      virtual ~ConstraintInterface() {}

      //! Copy
      virtual void copy(const ConstraintInterface& source) = 0;

      //! Cloning function
      virtual 
      Teuchos::RCP<LOCA::MultiContinuation::ConstraintInterface>
      clone(NOX::CopyType type = NOX::DeepCopy) const = 0;

      //! Return number of constraints
      virtual int numConstraints() const = 0;

      //! Set the solution vector to x.
      virtual void setX(const NOX::Abstract::Vector& x) = 0;

      //! Sets parameter indexed by paramID
      virtual void setParam(int paramID, double val) = 0;

      //! Sets parameters indexed by paramIDs
      virtual void setParams(const std::vector<int>& paramIDs, 
		     const NOX::Abstract::MultiVector::DenseMatrix& vals) = 0;

      //! Compute constraint residuals
      virtual NOX::Abstract::Group::ReturnType
      computeConstraints() = 0;

      //! Compute derivative of constraints w.r.t. solution vector x
      virtual NOX::Abstract::Group::ReturnType
      computeDX() = 0;

      //! Compute derivative of constraints w.r.t. supplied parameters.
      /*!
       * The first column of \c dgdp should be filled with the constraint
       * residuals \f$g\f$ if \c isValidG is \c false.  If \c isValidG is
       * \c true, then the \c dgdp contains \f$g\f$ on input.
       */
      virtual NOX::Abstract::Group::ReturnType
      computeDP(const std::vector<int>& paramIDs, 
		NOX::Abstract::MultiVector::DenseMatrix& dgdp, 
		bool isValidG) = 0;

      //! Return \c true if constraint residuals are valid
      virtual bool isConstraints() const = 0;

      //! Return \c true if derivative of constraint w.r.t. x is valid
      virtual bool isDX() const = 0;

      //! Return constraint residuals
      virtual const NOX::Abstract::MultiVector::DenseMatrix&
      getConstraints() const = 0;

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
	      NOX::Abstract::MultiVector::DenseMatrix& result_p) const = 0;

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
	      NOX::Abstract::MultiVector& result_x) const = 0;

      /*! 
       * \brief Return \c true if solution component of constraint 
       * derivatives is zero
       */
      virtual bool isDXZero() const = 0;

      //! Perform any preprocessing before a continuation step starts.
      /*!
       * The \c stepStatus argument indicates whether the previous step was
       * successful.  The default implementation is empty.
       */
      virtual void 
      preProcessContinuationStep(
			   LOCA::Abstract::Iterator::StepStatus stepStatus) {}
	
      //! Perform any postprocessing after a continuation step finishes.
      /*!
       * The \c stepStatus argument indicates whether the step was
       * successful.  The default implementation is empty.
       */
      virtual void 
      postProcessContinuationStep(
			   LOCA::Abstract::Iterator::StepStatus stepStatus) {}

    }; // Class ConstraintInterface

  } // Namespace MultiContinuation

} // Namespace LOCA

#endif // LOCA_MULTICONTINUATION_CONSTRAINTINTERFACE_H
