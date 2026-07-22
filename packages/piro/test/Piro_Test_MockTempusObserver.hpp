// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_TEST_MOCKTEMPUSOBSERVER_HPP
#define PIRO_TEST_MOCKTEMPUSOBSERVER_HPP

#include "Piro_Test_MockObserver.hpp"
#include "Tempus_IntegratorObserverBasic.hpp"

namespace Piro {

namespace Test {

template <typename Scalar>
class MockTempusObserver : public MockObserver<Scalar>,
                           public Tempus::IntegratorObserverBasic<Scalar>
{
public:
  MockTempusObserver() = default;

  void observeEndTimeStep(const Tempus::Integrator<Scalar>& integrator) override;
};

// ------------ IMPL ------------ //

template <typename Scalar>
void MockTempusObserver<Scalar>::
observeEndTimeStep(const Tempus::Integrator<Scalar>& integrator)
{
  // Simply grab stuff from integrator and forward to the base class observeSolution method
  // NOTE: this simple class IGNORES completely dxdp

  auto sol_hist = integrator.getSolutionHistory();

  auto x = sol_hist->getCurrentState()->getX();
  auto x_dot = sol_hist->getCurrentState()->getXDot();
  auto x_dotdot = sol_hist->getCurrentState()->getXDotDot();

  const Scalar scalar_time = sol_hist->getCurrentState()->getTime();
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType StampScalar;
  const StampScalar time = Teuchos::ScalarTraits<Scalar>::real(scalar_time);

  if (x_dot != Teuchos::null) {
    if (x_dotdot != Teuchos::null) {
      // x_dot AND x_dotdot
      this->observeSolution(*x, *x_dot, *x_dotdot, time);
    } else {
      // no x_dotdot
      this->observeSolution(*x, *x_dot, time);
    }
  } else {
    //no x_dot
    this->observeSolution(*x, time);
  }
}

} // namespace Test

} // namespace Piro

#endif /* PIRO_TEST_MOCKTEMPUSOBSERVER_HPP */

