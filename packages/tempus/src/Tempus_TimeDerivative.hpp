// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_TimeDerivative_hpp
#define Tempus_TimeDerivative_hpp

// Thrya
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"


namespace Tempus {

/** \brief This interface defines the time derivative connection between
 *         an implicit Stepper and WrapperModelEvaluator.
 *
 *  For a WrapperModelEvaluator which uses an implicit Stepper, i.e.,
 *  uses the implicit ODE or DAE form, \f$\mathcal{F}(\ddot{x},\dot{x},x,t)\f$,
 *  requires definition of the time derivatives during the iterations of
 *  the nonlinear solve. Note that if the Stepper solves for a time
 *  derivative, e.g., \f$\ddot{x}\f$ for the Newmark-\f$\beta\f$ methods,
 *  definitions for \f$x\f$ and \f$\dot{x} are required in the function
 *  compute().  This interface defines the calling function to compute
 *  those derivatives and/or state.
 */
template <typename Scalar>
class TimeDerivative
{
public:

  /// Set the underlying application ModelEvaluator
  virtual void compute(
    Teuchos::RCP<const Thyra::VectorBase<Scalar> > x,
    Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDot,
    Teuchos::RCP<      Thyra::VectorBase<Scalar> > xDotDot = Teuchos::null) = 0;

  // Derived classes may need an initialize, but the argument lists will vary.
  // virtual void initialize(Scalar dt, ... ) = 0;
};


} // namespace Tempus
#endif // Tempus_TimeDerivative_hpp
