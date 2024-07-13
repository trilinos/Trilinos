// @HEADER
// *****************************************************************************
//            LOCA: Library of Continuation Algorithms Package
//
// Copyright 2001-2005 NTESS and the LOCA contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef LOCA_STEPSIZE_ADAPTIVE_H
#define LOCA_STEPSIZE_ADAPTIVE_H

#include "LOCA_StepSize_Constant.H"  // base class

namespace LOCA {

  namespace StepSize {

    //! %Adaptive step size control strategy
    /*!
     * This class implements an adaptive step size control strategy derived
     * from the strategy implemented in the LOCA::StepSize::Constant class.
     * If the previous step was unsucessful, the step size is cut in half as
     * in the constant strategy, but if the step was sucessful this strategy
     * increases the step size based on the number of nonlinear solver
     * iterations required in the previous step.  In particular, the new
     * step size \f$\Delta s_{new}\f$ is given by
     * \f[
     *    \Delta s_{new} = \Delta s_{old}\left(1 + a\left(\frac{N_{max} - N}{N_{max}}\right)^2\right)
     * \f]
     * where \f$a\in[0,1]\f$ is an aggressiveness factor, \f$N\f$ is the
     * number of nonlinear solver iterations in the previous step, and
     * \f$N_{max}\f$ is the maximum number of nonlinear solver iterations.
     *
     * The parameters used by this class supplied in the constructor are
     * the same as used by the Constant class in addition to:
     * <ul>
     * <li> "Aggressiveness" - Aggressiveness factor \f$a\f$ (Default 0.5)
     * </ul>
     */
    class Adaptive : public LOCA::StepSize::Constant {

    public:

      //! Constructor.
      Adaptive(
     const Teuchos::RCP<LOCA::GlobalData>& global_data,
     const Teuchos::RCP<LOCA::Parameter::SublistParser>& topParams,
     const Teuchos::RCP<Teuchos::ParameterList>& stepsizeParams);

      //! Destructor
      virtual ~Adaptive();

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
//              const LOCA::Stepper& stepper,
              const LOCA::Abstract::Iterator& stepper,
              double& stepSize);

    protected:

      //! Stores the aggressiveness factor \f$a\f$
      double agrValue;

      //! Maximum number of nonlinear iterations
      double maxNonlinearSteps;
    };
  } // namespace StepSize
} // namespace LOCA

#endif
