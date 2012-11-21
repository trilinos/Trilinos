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

#ifndef LOCA_STEPSIZE_CONSTANT_H
#define LOCA_STEPSIZE_CONSTANT_H

#include "Teuchos_RCP.hpp" 

#include "LOCA_StepSize_AbstractStrategy.H"  // base class

// forward declarations
namespace LOCA {
  class GlobalData;
  namespace Parameter {
    class SublistParser;
  }
}

namespace LOCA {

  namespace StepSize {

    //! %Constant step size control strategy
    /*!
     * This class implements a roughly constant step size control strategy.
     * If the previous step was sucessful, the new step size is set equal
     * to the old, otherwise the step size is cut by a supplied factor.  
     * Once a sucessful
     * step is made, the step size is increased by a supplied factor 
     * until the initial step size is reached.
     *
     * This class also incorporates rescaling of the continuation parameter
     * when calculating a step size (common in arc-length continuation).
     * For the first continuation step, the step size is chosen so that
     * step size times the parameter component of the predictor is equal
     * to the initial step size.  From then on, the step size is multiplied
     * by the step size scale factor (see 
     ( LOCA::MultiContinuation::ArcLengthGroup)
     * which incorporates rescaling of the continuation parameter.
     *
     * The parameters used by this class supplied in the constructor are:
     * <ul>
     * <li> "Max Step Size" - Largest valid step size (Default 1.0e+12)
     * <li> "Min Step Size" - Smallest valid step size (Default 1.0e-12)
     * <li> "Initial Step Size" - Initial step size (Default 1.0)
     * <li> "Failed Step Reduction Factor" - Factor by which step size is 
     *      reduced after a failed step. (Default 0.5)
     * <li> "Successful Step Increase Factor" - Factor by which step size is
     *      increased after a successful step when the step size is smaller
     *      than the initial step size. (Default 1.26 = \f$2^{1/3}\f$)
     * </ul>
     */
    class Constant : public LOCA::StepSize::AbstractStrategy {

    public:

      //! Constructor. 
      Constant(
         const Teuchos::RCP<LOCA::GlobalData>& global_data,
	 const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
	 const Teuchos::RCP<Teuchos::ParameterList>& stepsizeParams);

      //! Destructor
      virtual ~Constant();

      //! Compute the step size as described above
      /*! 
       * \param curGroup [in] Current continuation group
       * \param predictor [in] Current predictor direction
       * \param solver [in] Solver from previous step
       * \param stepStatus [in] Status of previous step
       * \param stepper [in] Stepper
       * \param stepSize [out] Computed step size
       * \returns ReturnType code indicating success or failure
       * Returns NOX::Abstract::Group::Failed if the computed step size is
       * smaller than the minimum step size
       */
      virtual NOX::Abstract::Group::ReturnType 
      computeStepSize(LOCA::MultiContinuation::AbstractStrategy& curGroup,
		      const LOCA::MultiContinuation::ExtendedVector& predictor,
		      const NOX::Solver::Generic& solver,
		      const LOCA::Abstract::Iterator::StepStatus& stepStatus,
//		      const LOCA::Stepper& stepper,
		      const LOCA::Abstract::Iterator& stepper,
		      double& stepSize);

      //! Returns previous step size
      virtual double getPrevStepSize() const;

      //! Returns initial step size
      virtual double getStartStepSize() const;

    protected:

      /*! 
       * \brief Clip the computed step size to the bounds given by the 
       * maximum and minimum step sizes
       */
      /*! 
       * Returns NOX::Abstract::Group::Failed if the computed step size is
       * smaller than the minimum step size
       */
      virtual NOX::Abstract::Group::ReturnType 
      clipStepSize(double& stepSize);

    protected:

      //! Global data
      Teuchos::RCP<LOCA::GlobalData> globalData;

      //! Maximum step size
      double maxStepSize;

      //! Minimum step size
      double minStepSize;

      //! Initial step size
      double startStepSize;

      //! Factor by which step size is reduced after a failed step
      double failedFactor;

      //! Factor by which step size is increased after a successful step
      double successFactor;

      //! Previous step size
      double prevStepSize;

      //! Flag indicating if this is the first step
      bool isFirstStep;
    };
  } // namespace StepSize
} // namespace LOCA

#endif
