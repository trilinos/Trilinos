// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_WrapperModelEvaluator_hpp
#define Tempus_WrapperModelEvaluator_hpp

#include <functional>
#include "Thyra_StateFuncModelEvaluatorBase.hpp"

namespace Tempus {

/** \brief A ModelEvaluator which wraps the application ModelEvaluator.
 *
 *  The WrapperModelEvaluator takes a state, \f$x\f$, computes time
 *  derivative(s), \f$\dot{x}\f$ and/or \f$\ddot{x}\f$, from the
 *  implicit stepper (StepperImplicit) and calls the application
 *  ModelEvaluator to determine its residual, \f$\mathcal{F}(x)\f$,
 *  which is suitable for the nonlinear solve.
 */
template <typename Scalar>
class WrapperModelEvaluator : public Thyra::StateFuncModelEvaluatorBase<Scalar>
{
public:

  /// Set the underlying application ModelEvaluator
  virtual void setAppModel(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > & me) = 0;

  /// Get the underlying application model 'f'
  virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
    getAppModel() const = 0;

  /// Set values to compute x dot and evaluate application model.
  virtual void initialize(
    std::function<void (const Thyra::VectorBase<Scalar> &,
                              Thyra::VectorBase<Scalar> &)> computeXDot,
    double t, double alpha, double beta) = 0;
};

} // namespace Tempus

#endif // Tempus_WrapperModelEvaluator_hpp
