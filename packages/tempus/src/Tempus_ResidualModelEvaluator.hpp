// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_ResidualModelEvaluator_hpp
#define Tempus_ResidualModelEvaluator_hpp

#include <functional>
#include "Thyra_StateFuncModelEvaluatorBase.hpp"

namespace Tempus {

/** \brief A ModelEvaluator for residual evaluations given a state.
 *  This ModelEvaluator takes a state, x, and determines its residual,
 *  \f$ g(x) \f$, which is suitable for a nonlinear solve.  This is
 *  accomplished by computing the time derivative of the state, x_dot,
 *  (through Lambda functions), supplying the current time, and calling
 *  the application transient ModelEvaluator, \f$ f(\dot{x},x,t) \f$.
 *
 *  This class breaks the primary design principle for ModelEvaluators;
 *  it is not stateless!
 */
template <typename Scalar>
class ResidualModelEvaluator : public Thyra::StateFuncModelEvaluatorBase<Scalar>
{
public:

  /// Set the underlying transient ModelEvaluator
  virtual void setTransientModel(
    const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> > & me) = 0;

  /// Get the underlying transient model 'f'
  virtual Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >
    getTransientModel() const = 0;

  /// Set values to compute x dot and evaluate transient model.
  virtual void initialize(
    std::function<void (const Thyra::VectorBase<Scalar> &,
                              Thyra::VectorBase<Scalar> &)> computeXDot,
    double t, double alpha, double beta) = 0;
};

} // namespace Tempus

#endif // Tempus_ResidualModelEvaluator_hpp
