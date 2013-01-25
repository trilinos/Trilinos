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

#ifndef LOCA_MULTIPREDICTOR_ABSTRACTSTRATEGY_H
#define LOCA_MULTIPREDICTOR_ABSTRACTSTRATEGY_H

#include "Teuchos_RCP.hpp"

#include "NOX_Abstract_Group.H"  // for ReturnType, Teuchos::ParameterList

// Forward declarations
namespace LOCA {
  namespace MultiContinuation {
    class ExtendedGroup;
    class ExtendedVector;
    class ExtendedMultiVector;
  }
}

namespace LOCA {

  //! Predictor direction strategies
  namespace MultiPredictor {
    
    //! Abstract interface class for predictor strategies
    /*!
     * AbstractStrategy defines an abstract interface for predictor strategies.
     * It is used by the LOCA::Stepper and LOCA::MultiContinuation groups to
     * compute predictors and tangent approximations to continuation curves
     * and surfaces.
     *
     * The interface defines several pure virtual methods that derived classes
     * should implement for a particular predictor strategy.  Note that
     * predictor strategies are assumed to have a state, and therefore need
     * to define copy constructors, assignment operators, and a clone
     * function.  Constructors of derived classes should be of the form:
     *
     * \code
     * class Derived : public AbstractStrategy {
     * public:
     *   Derived(
     *   const Teuchos::RCP<LOCA::GlobalData>& global_data,
     *   const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& predictorParams);
     *   ...
     * };
     * \endcode
     *
     * where \c global_data is the LOCA global data object, \c topParams is
     * the parsed top-level parameter list, and \c predictorParams 
     * is a parameter list of predictor parameters.
     *
     * This class and its children follow the Strategy pattern as defined
     * in Erich Gamma, et al. "Design Patterns:  Elements of Reusable 
     * Object-Oriented Software." Addison Wesley, Boston, MA, 1995.
     */
    class AbstractStrategy {

    public:

      //! Constructor. 
      AbstractStrategy() {}

      //! Destructor
      virtual ~AbstractStrategy() {}

      //! Assignment operator
      virtual LOCA::MultiPredictor::AbstractStrategy& 
      operator=(const LOCA::MultiPredictor::AbstractStrategy& source) = 0;

      //! Clone function
      virtual Teuchos::RCP<LOCA::MultiPredictor::AbstractStrategy>
      clone(NOX::CopyType type = NOX::DeepCopy) const = 0;

      /*! 
       * \brief Compute the predictor given the current and
       * previous solution vectors.  Set \em baseOnSecant to \em false
       * if the predictor orientation should not be based on the
       * secant vector (first or last steps of a continuation run).
       */
      /*!
       * As an example for a first-order predictor, this method should
       * compute the approximate tangent to the continuation curve.
       */
      virtual NOX::Abstract::Group::ReturnType 
      compute(bool baseOnSecant, const std::vector<double>& stepSize,
	      LOCA::MultiContinuation::ExtendedGroup& grp,
	      const LOCA::MultiContinuation::ExtendedVector& prevXVec,
	      const LOCA::MultiContinuation::ExtendedVector& xVec) = 0;

      //! Evaluate predictor with step size \c stepSize.
      /*!
       * For a first-order predictor, this method should compute
       * \c result[i] = \c xVec[i] + \c stepSize[i] * v[i]
       * for each i, where \c v[i] is the ith predictor direction.
       */
      virtual NOX::Abstract::Group::ReturnType 
      evaluate(const std::vector<double>& stepSize,
	       const LOCA::MultiContinuation::ExtendedVector& xVec,
	       LOCA::MultiContinuation::ExtendedMultiVector& result) const = 0;

      //! Compute tangent to predictor and store in \c tangent.
      /*!
       * For a first-order predictor, this is the predictor direction itself.
       */
      virtual NOX::Abstract::Group::ReturnType 
      computeTangent(LOCA::MultiContinuation::ExtendedMultiVector& tangent)= 0;

      //! Is the tangent vector for this predictor scalable
      /*!
       * This method determines whether the approximate tangent computed
       * by this strategy is appropriate for scaling.
       */
      virtual bool isTangentScalable() const = 0;

    protected:

      /*! 
       * \brief Sets orientation of predictor based on parameter change 
       * from previous steps
       */
      /*!
       * The implementation here looks at the sign of the scaled dot product
       * between the secant vector and the predictor and changes the sign
       * of the predictor if this scaled dot product is a different sign
       * than the step size.  If \em baseOnSecant is \em false, then the
       * sign is chosen so the parameter component of the predictor is 
       * positive for cases when a secant vector is not available (first
       * step in a continuation run) or may give incorrect information
       * (last step of a continuation run).
       */
      virtual void
      setPredictorOrientation(
	             bool baseOnSecant, 
		     const std::vector<double>& stepSize,
		     const LOCA::MultiContinuation::ExtendedGroup& grp,
	             const LOCA::MultiContinuation::ExtendedVector& prevXVec,
		     const LOCA::MultiContinuation::ExtendedVector& xVec,
		     LOCA::MultiContinuation::ExtendedVector& secant,
	             LOCA::MultiContinuation::ExtendedMultiVector& tangent);

    }; // Class AbstractStrategy

  } // Namespace MultiPredictor

} // Namespace LOCA

#endif
