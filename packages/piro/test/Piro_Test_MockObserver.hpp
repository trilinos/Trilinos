// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_TEST_MOCKOBSERVER_HPP
#define PIRO_TEST_MOCKOBSERVER_HPP

#include "Piro_ObserverBase.hpp"

#include "Thyra_VectorBase.hpp"

namespace Piro {

namespace Test {

template <typename Scalar>
class MockObserver : public ObserverBase<Scalar> {
public:
  MockObserver();

  virtual void observeSolution(
      const Thyra::VectorBase<Scalar> &solution);

  virtual void observeSolution(
      const Thyra::VectorBase<Scalar> &solution,
      const Thyra::MultiVectorBase<Scalar> &solution_dxdp); 
  
  virtual void observeSolution(
      const Thyra::VectorBase<Scalar> &solution,
      const Scalar stamp);
  
  virtual void observeSolution(
      const Thyra::VectorBase<Scalar> &solution,
      const Thyra::MultiVectorBase<Scalar> &solution_dxdp,
      const Scalar stamp);

  virtual void observeSolution(
      const Thyra::VectorBase<Scalar> &solution,
      const Thyra::VectorBase<Scalar> &solution_dot,
      const Scalar stamp);
  
  virtual void observeSolution(
      const Thyra::VectorBase<Scalar> &solution,
      const Thyra::MultiVectorBase<Scalar> &solution_dxdp,
      const Thyra::VectorBase<Scalar> &solution_dot,
      const Scalar stamp);
  
  virtual void observeSolution(
      const Thyra::VectorBase<Scalar> &solution,
      const Thyra::VectorBase<Scalar> &solution_dot,
      const Thyra::VectorBase<Scalar> &solution_dotdot,
      const Scalar stamp);
  
  virtual void observeSolution(
      const Thyra::VectorBase<Scalar> &solution,
      const Thyra::MultiVectorBase<Scalar> &solution_dxdp,
      const Thyra::VectorBase<Scalar> &solution_dot,
      const Thyra::VectorBase<Scalar> &solution_dotdot,
      const Scalar stamp);

  Teuchos::RCP<const Thyra::VectorBase<Scalar> > lastSolution() const {
    return lastSolution_;
  }

  Teuchos::RCP<const Thyra::MultiVectorBase<Scalar> > lastSolution_dxdp() const {
    return lastSolution_dxdp_;
  }

  Teuchos::RCP<const Thyra::VectorBase<Scalar> > lastSolution_dot() const {
    return lastSolution_dot_;
  }

  Teuchos::RCP<const Thyra::VectorBase<Scalar> > lastSolution_dotdot() const {
    return lastSolution_dotdot_;
  }

  Scalar lastStamp() const {
    return lastStamp_;
  }

private:
  Teuchos::RCP<Thyra::VectorBase<Scalar> > lastSolution_;
  Teuchos::RCP<Thyra::MultiVectorBase<Scalar> > lastSolution_dxdp_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > lastSolution_dot_;
  Teuchos::RCP<Thyra::VectorBase<Scalar> > lastSolution_dotdot_;
  Scalar lastStamp_;
};


template <typename Scalar>
MockObserver<Scalar>::MockObserver() :
  lastSolution_(Teuchos::null),
  lastSolution_dxdp_(Teuchos::null),
  lastSolution_dot_(Teuchos::null),
  lastSolution_dotdot_(Teuchos::null),
  lastStamp_(Scalar())
{}

template <typename Scalar>
void
MockObserver<Scalar>::observeSolution(
    const Thyra::VectorBase<Scalar> &solution)
{
  lastSolution_ = solution.clone_v();
}

template <typename Scalar>
void
MockObserver<Scalar>::observeSolution(
    const Thyra::VectorBase<Scalar> &solution,
    const Thyra::MultiVectorBase<Scalar> &solution_dxdp)
{
  lastSolution_ = solution.clone_v();
  lastSolution_dxdp_ = solution_dxdp.clone_mv();
}


template <typename Scalar>
void
MockObserver<Scalar>::observeSolution(
    const Thyra::VectorBase<Scalar> &solution,
    const Scalar stamp)
{
  lastSolution_ = solution.clone_v();
  lastStamp_ = stamp;
}

template <typename Scalar>
void
MockObserver<Scalar>::observeSolution(
    const Thyra::VectorBase<Scalar> &solution,
    const Thyra::MultiVectorBase<Scalar> &solution_dxdp,
    const Scalar stamp)
{
  lastSolution_ = solution.clone_v();
  lastSolution_dxdp_ = solution_dxdp.clone_mv();
  lastStamp_ = stamp;
}

template <typename Scalar>
void
MockObserver<Scalar>::observeSolution(
    const Thyra::VectorBase<Scalar> &solution,
    const Thyra::VectorBase<Scalar> &solution_dot,
    const Scalar stamp)
{
  lastSolution_ = solution.clone_v();
  lastSolution_dot_ = solution_dot.clone_v();
  lastStamp_ = stamp;
}

template <typename Scalar>
void
MockObserver<Scalar>::observeSolution(
    const Thyra::VectorBase<Scalar> &solution,
    const Thyra::MultiVectorBase<Scalar> &solution_dxdp,
    const Thyra::VectorBase<Scalar> &solution_dot,
    const Scalar stamp)
{
  lastSolution_ = solution.clone_v();
  lastSolution_dxdp_ = solution_dxdp.clone_mv();
  lastSolution_dot_ = solution_dot.clone_v();
  lastStamp_ = stamp;
}

template <typename Scalar>
void
MockObserver<Scalar>::observeSolution(
    const Thyra::VectorBase<Scalar> &solution,
    const Thyra::VectorBase<Scalar> &solution_dot,
    const Thyra::VectorBase<Scalar> &solution_dotdot,
    const Scalar stamp)
{
  lastSolution_ = solution.clone_v();
  lastSolution_dot_ = solution_dot.clone_v();
  lastSolution_dotdot_ = solution_dotdot.clone_v();
  lastStamp_ = stamp;
}

template <typename Scalar>
void
MockObserver<Scalar>::observeSolution(
    const Thyra::VectorBase<Scalar> &solution,
    const Thyra::MultiVectorBase<Scalar> &solution_dxdp, 
    const Thyra::VectorBase<Scalar> &solution_dot,
    const Thyra::VectorBase<Scalar> &solution_dotdot,
    const Scalar stamp)
{
  lastSolution_ = solution.clone_v();
  lastSolution_dxdp_ = solution_dxdp.clone_mv();
  lastSolution_dot_ = solution_dot.clone_v();
  lastSolution_dotdot_ = solution_dotdot.clone_v();
  lastStamp_ = stamp;
}

} // namespace Test

} // namespace Piro

#endif /* PIRO_TEST_MOCKOBSERVER_HPP */
