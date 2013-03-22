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

#ifndef LOCA_STEPSIZE_ABSTRACTSTRATEGY_H
#define LOCA_STEPSIZE_ABSTRACTSTRATEGY_H

#include "NOX_Abstract_Group.H"     // for ReturnType
#include "LOCA_Abstract_Iterator.H" // for StepStatus
#include "Teuchos_ParameterList.hpp"

// Forward declarations
namespace NOX {
  namespace Solver {
    class Generic;
  }
}
namespace LOCA {
  namespace MultiContinuation {
    class AbstractStrategy;
    class ExtendedVector;
  }
//  class Stepper;
}

namespace LOCA {

  //! Collection of step size control strategies
  namespace StepSize {
    
    //! Abstract interface class for step size control strategies
    /*!
     * AbstractStrategy defines an abstract interface for step size control 
     * strategies. It is used by LOCA::Stepper to the step size for each 
     * continuation step.  
     *
     * The interface currently defines three pure virtual methods, 
     * computeStepSize() to compute the step size, getPrevSteSize() to get
     * the step size from the previous step, and getStartStepSize() to get 
     * the initial step size.  Derived classes should implement this method 
     * for a particular strategy.  Constructors for derived classes should be 
     * of the form:
     *
     * \code
     * class Derived : public AbstractStrategy {
     * public:
     *   Derived(
     *   const Teuchos::RCP<LOCA::GlobalData>& global_data,
     *   const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& stepsizeParams);
     *   ...
     * };
     * \endcode
     *
     * where \c global_data is the LOCA global data object, \c topParams is
     * the parsed top-level parameter list, and \c stepsizeParams 
     * is a parameter list of step size control parameters.
     *
     * This class and its children follow the Strategy pattern as defined
     * in Erich Gamma, et al. "Design Patterns:  Elements of Reusable 
     * Object-Oriented Software." Addison Wesley, Boston, MA, 1995.
     */
    class AbstractStrategy {

    public:

      //! Constructor
      AbstractStrategy() {}

      //! Destructor
      virtual ~AbstractStrategy() {}

      //! Compute step size
      /*!
       * \param curGroup [in] Current continuation group
       * \param predictor [in] Current predictor direction
       * \param solver [in] Solver from previous step
       * \param stepStatus [in] Status of previous step
       * \param stepper [in] Stepper
       * \param stepSize [out] Computed step size
       * \returns ReturnType code indicating success or failure
       */
      virtual NOX::Abstract::Group::ReturnType 
      computeStepSize(LOCA::MultiContinuation::AbstractStrategy& curGroup,
		      const LOCA::MultiContinuation::ExtendedVector& predictor,
		      const NOX::Solver::Generic& solver,
		      const LOCA::Abstract::Iterator::StepStatus& stepStatus,
//		      const LOCA::Stepper& stepper,
		      const LOCA::Abstract::Iterator& stepper,
		      double& stepSize) = 0;

      //! Return the previous step size
      virtual double getPrevStepSize() const = 0;

      //! Return the initial step size
      virtual double getStartStepSize() const = 0;

    private:

      //! Private to prohibit copying
      AbstractStrategy(const AbstractStrategy&);

      //! Private to prohibit copying
      AbstractStrategy& operator = (const AbstractStrategy&);

    }; // Class AbstractStrategy

  } // Namespace StepSize

} // Namespace LOCA

#endif // LOCA_STEPSIZE_ABSTRACTSTRATEGY_H
