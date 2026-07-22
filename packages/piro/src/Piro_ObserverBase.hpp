// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIRO_OBSERVERBASE_HPP
#define PIRO_OBSERVERBASE_HPP

#include <string>

#include "Thyra_VectorBase.hpp"
#include "Thyra_DefaultMultiVectorProductVector.hpp"

namespace Piro {

template <typename Scalar>
class ObserverBase {
public:
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

  virtual void observeSolution(
      const Thyra::MultiVectorBase<Scalar> &solution, Scalar time);

  virtual void observeSolution(
      const Thyra::MultiVectorBase<Scalar> &solution, 
      const Thyra::MultiVectorBase<Scalar> &solution_dxdp,
      Scalar time);

  virtual void parameterChanged(
      const std::string& param);

  virtual ~ObserverBase() {}
};

template <typename Scalar>
void
ObserverBase<Scalar>::observeSolution(
    const Thyra::VectorBase<Scalar> &/*solution*/)
{
  // Nothing to do by default
}

template <typename Scalar>
void
ObserverBase<Scalar>::observeSolution(
    const Thyra::VectorBase<Scalar> &/*solution*/,
    const Thyra::MultiVectorBase<Scalar> &/*solution_dxdp*/)
{
  // Nothing to do by default
}

template <typename Scalar>
void
ObserverBase<Scalar>::observeSolution(
    const Thyra::VectorBase<Scalar> &/*solution*/,
    const Scalar /*stamp*/)
{
  // Nothing to do by default
}

template <typename Scalar>
void
ObserverBase<Scalar>::observeSolution(
    const Thyra::VectorBase<Scalar> &/*solution*/,
    const Thyra::MultiVectorBase<Scalar> &/*solution_dxdp*/,
    const Scalar /*stamp*/)
{
  // Nothing to do by default
}

template <typename Scalar>
void
ObserverBase<Scalar>::observeSolution(
    const Thyra::VectorBase<Scalar> &/*solution*/,
    const Thyra::VectorBase<Scalar> &/*solution_dot*/,
    const Scalar /*stamp*/)
{
  // Nothing to do by default
}

template <typename Scalar>
void
ObserverBase<Scalar>::observeSolution(
    const Thyra::VectorBase<Scalar> &/*solution*/,
    const Thyra::MultiVectorBase<Scalar> &/*solution_dxdp*/,
    const Thyra::VectorBase<Scalar> &/*solution_dot*/,
    const Scalar /*stamp*/)
{
  // Nothing to do by default
}

template <typename Scalar>
void
ObserverBase<Scalar>::observeSolution(
    const Thyra::VectorBase<Scalar> &/*solution*/,
    const Thyra::VectorBase<Scalar> &/*solution_dot*/,
    const Thyra::VectorBase<Scalar> &/*solution_dotdot*/,
    const Scalar /*stamp*/)
{
  // Nothing to do by default
}

template <typename Scalar>
void
ObserverBase<Scalar>::observeSolution(
    const Thyra::VectorBase<Scalar> &/*solution*/,
    const Thyra::MultiVectorBase<Scalar> &/*solution_dxdp*/,
    const Thyra::VectorBase<Scalar> &/*solution_dot*/,
    const Thyra::VectorBase<Scalar> &/*solution_dotdot*/,
    const Scalar /*stamp*/)
{
  // Nothing to do by default
}

template <typename Scalar>
void
ObserverBase<Scalar>::observeSolution(
      const Thyra::MultiVectorBase<Scalar> &/*solution*/, 
      Scalar /*time*/)
{
  // Nothing to do by default
}

template <typename Scalar>
void
ObserverBase<Scalar>::observeSolution(
      const Thyra::MultiVectorBase<Scalar> &/*solution*/, 
      const Thyra::MultiVectorBase<Scalar> &/*solution_dxdp*/,
      Scalar /*time*/)
{
  // Nothing to do by default
}

template <typename Scalar>
void
ObserverBase<Scalar>::parameterChanged(
      const std::string& param)
{
  // Nothing to do by default
}

} // namespace Piro

#endif /* PIRO_OBSERVERBASE_HPP */
