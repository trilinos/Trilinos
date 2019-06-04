// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_SensitivityModelEvaluatorBase_hpp
#define Tempus_SensitivityModelEvaluatorBase_hpp

#include "Thyra_ModelEvaluatorDefaultBase.hpp"
#include "Thyra_NonlinearSolverBase.hpp"
#include "Tempus_SolutionHistory.hpp"

namespace Tempus {

/** \brief A ModelEvaluator decorator for sensitivity analysis
 *
 * Used in sensitivity analysis model evaluators for interpolating from a
 * previous solution and such.  All additional methods have default, empty
 * implementations, since different sensitivity model evaluators need slightly
 * different capabilities.
 */
template <typename Scalar>
class SensitivityModelEvaluatorBase :
    public virtual Thyra::ModelEvaluatorDefaultBase<Scalar> {
public:

  //! Constructor
  SensitivityModelEvaluatorBase() {}

  //! Destructor
  virtual ~SensitivityModelEvaluatorBase() {}

  //! Get the underlying forward model
  virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
  getForwardModel() const { return Teuchos::null; }

  //! Set solution history from forward state evaluation (for interpolation)
  virtual void setForwardSolutionHistory(
    const Teuchos::RCP<const Tempus::SolutionHistory<Scalar> >& /* sh */) {}

  //! Set solution state from forward state evaluation (for frozen state)
    virtual void setForwardSolutionState(
      const Teuchos::RCP<const Tempus::SolutionState<Scalar> >& /* s */) {}

  //! Set the solver of the underlying model if you want to reuse it
  virtual void setSolver(
    const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& /* solver */,
    const bool /* force_W_update */) {}
};

} // namespace Tempus

#endif
