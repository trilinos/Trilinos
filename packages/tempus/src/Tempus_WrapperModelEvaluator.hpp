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

/** \brief A ModelEvaluator for residual evaluations given a state.
 *  This ModelEvaluator takes a state, x, and determines its residual,
 *  \f$ g(x) \f$, which is suitable for a nonlinear solve.  This is
 *  accomplished by computing the time derivative of the state, x_dot,
 *  (through Lambda functions), supplying the current time, and calling
 *  the application ModelEvaluator, \f$ f(\dot{x},x,t) \f$.
 *
 *  This class breaks the primary design principle for ModelEvaluators;
 *  it is not stateless!
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
