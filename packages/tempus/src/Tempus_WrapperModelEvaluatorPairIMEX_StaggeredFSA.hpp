//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_ModelEvaluatorPairIMEX_StaggeredFSA_hpp
#define Tempus_ModelEvaluatorPairIMEX_StaggeredFSA_hpp

#include "Tempus_config.hpp"
#include "Tempus_SensitivityModelEvaluatorBase.hpp"
#include "Tempus_WrapperModelEvaluatorPairIMEX_Basic.hpp"
#include "Tempus_StaggeredForwardSensitivityModelEvaluator.hpp"

namespace Tempus {

/** \brief Specialization of IMEX ME for "staggered" FSA method.
 *
 * This specializes the implementation of several parts of
 * WrapperModelEvaluatorPairIMEX_Basic for forward-sensitivity analysis
 * with StaggeredForwardSensitivityModelEvaluator.
 */
template <typename Scalar>
class WrapperModelEvaluatorPairIMEX_StaggeredFSA
  : public SensitivityModelEvaluatorBase<Scalar>,
    public WrapperModelEvaluatorPairIMEX_Basic<Scalar> {
 public:
  /// Constructor
  WrapperModelEvaluatorPairIMEX_StaggeredFSA(
      const Teuchos::RCP<const WrapperModelEvaluatorPairIMEX_Basic<Scalar> >&
          forwardModel,
      const bool is_pseudotransient,
      const Teuchos::RCP<const Teuchos::ParameterList>& pList = Teuchos::null)
  {
    forwardModel_     = forwardModel;
    appExplicitModel_ = forwardModel_->getExplicitModel();
    appImplicitModel_ = forwardModel_->getImplicitModel();
    fsaExplicitModel_ =
        rcp(new FSAME(appExplicitModel_, appExplicitModel_, appExplicitModel_,
                      is_pseudotransient, pList));
    fsaImplicitModel_ =
        rcp(new FSAME(appImplicitModel_, appImplicitModel_, appImplicitModel_,
                      is_pseudotransient, pList));
    Base::setup(fsaExplicitModel_, fsaImplicitModel_);
  }

  /// Destructor
  virtual ~WrapperModelEvaluatorPairIMEX_StaggeredFSA() {}

  /// \name Overridden from Tempus::SensitivityModelEvaluatorBase
  //@{

  /// Get the underlying forward model
  virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getForwardModel()
      const
  {
    return forwardModel_;
  }

  /// Set solution history from forward state evaluation (for interpolation)
  virtual void setForwardSolutionHistory(
      const Teuchos::RCP<const Tempus::SolutionHistory<Scalar> >& sh)
  {
    fsaExplicitModel_->setForwardSolutionHistory(sh);
    fsaImplicitModel_->setForwardSolutionHistory(sh);
  }

  /// Set solution state from forward state evaluation (for frozen state)
  virtual void setForwardSolutionState(
      const Teuchos::RCP<const Tempus::SolutionState<Scalar> >& s)
  {
    fsaExplicitModel_->setForwardSolutionState(s);
    fsaImplicitModel_->setForwardSolutionState(s);
  }

  /// Set the solver of the underlying model if you want to reuse it
  virtual void setSolver(
      const Teuchos::RCP<Thyra::NonlinearSolverBase<Scalar> >& solver,
      const bool force_W_update)
  {
    fsaImplicitModel_->setSolver(solver, force_W_update);
  }

  //@}

 private:
  /// Default constructor - not allowed
  WrapperModelEvaluatorPairIMEX_StaggeredFSA() {}

 protected:
  typedef WrapperModelEvaluatorPairIMEX_Basic<Scalar> Base;
  typedef StaggeredForwardSensitivityModelEvaluator<Scalar> FSAME;

  Teuchos::RCP<const WrapperModelEvaluatorPairIMEX_Basic<Scalar> >
      forwardModel_;
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > appExplicitModel_;
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > appImplicitModel_;
  Teuchos::RCP<FSAME> fsaExplicitModel_;
  Teuchos::RCP<FSAME> fsaImplicitModel_;
};

}  // namespace Tempus

#endif  // Tempus_ModelEvaluatorPairIMEX_StaggeredFSA_hpp
