//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_AuxiliaryIntegralModelEvaluator_decl_hpp
#define Tempus_AuxiliaryIntegralModelEvaluator_decl_hpp

#include "Thyra_StateFuncModelEvaluatorBase.hpp"

#include "Tempus_config.hpp"
#include "Tempus_SolutionHistory.hpp"

namespace Tempus {

/** \brief ModelEvaluator for integrating auxiliary equations. */
/**
 * Given a ModelEvaluator defining an auxiliary/response function g(x_dot,x,p,t)
 * this ModelEvaluator defines a new ODE
 *     dz/dt = g(x_dot(t),x(t),p,t)
 * used to compute the integral of g over some time range.  This ModelEvaluator
 * can then be used by Tempus to compute the definite integral.
 */
template <typename Scalar>
class AuxiliaryIntegralModelEvaluator
  : public Thyra::StateFuncModelEvaluatorBase<Scalar> {
 public:
  typedef Thyra::VectorBase<Scalar> Vector;
  typedef Thyra::MultiVectorBase<Scalar> MultiVector;

  //! Constructor
  AuxiliaryIntegralModelEvaluator(
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > &model,
      const int g_index);

  //! Get the underlying model 'f'
  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > getModel() const
  {
    return model_;
  }

  //! Set solution history from forward evaluation
  void setForwardSolutionHistory(
      const Teuchos::RCP<const Tempus::SolutionHistory<Scalar> > &sh);

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_p_space(int p) const;

  Teuchos::RCP<const Teuchos::Array<std::string> > get_p_names(int p) const;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_x_space() const;

  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > get_f_space() const;

  Teuchos::RCP<Thyra::LinearOpBase<Scalar> > create_W_op() const;

  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar> >
  get_W_factory() const;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;

  //@}

 private:
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;

  void evalModelImpl(
      const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs,
      const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &outArgs) const;

  Thyra::ModelEvaluatorBase::InArgs<Scalar> prototypeInArgs_;
  Thyra::ModelEvaluatorBase::OutArgs<Scalar> prototypeOutArgs_;

  Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > model_;
  Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar> > space_;
  Teuchos::RCP<const Tempus::SolutionHistory<Scalar> > sh_;
  int g_index_;

  mutable Teuchos::RCP<Tempus::SolutionState<Scalar> > forward_state_;
  mutable Scalar t_interp_;
};

}  // namespace Tempus

#endif
