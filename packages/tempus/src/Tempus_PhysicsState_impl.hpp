//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhysicsState_impl_hpp
#define Tempus_PhysicsState_impl_hpp

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include <string>

namespace Tempus {

template <class Scalar>
PhysicsState<Scalar>::PhysicsState(std::string pN) : physicsName_(pN)
{
}

template <class Scalar>
Teuchos::RCP<PhysicsState<Scalar> > PhysicsState<Scalar>::clone() const
{
  Teuchos::RCP<PhysicsState<Scalar> > pSB =
      Teuchos::rcp(new PhysicsState<Scalar>());

  pSB->setName(this->getName());

  return pSB;
}

template <class Scalar>
void PhysicsState<Scalar>::copy(
    const Teuchos::RCP<const PhysicsState<Scalar> >& pS)
{
  physicsName_ = pS->getName();
}

template <class Scalar>
std::string PhysicsState<Scalar>::getName() const
{
  return physicsName_;
}

template <class Scalar>
void PhysicsState<Scalar>::setName(std::string pN)
{
  physicsName_ = pN;
}

template <class Scalar>
std::string PhysicsState<Scalar>::description() const
{
  return "Tempus::PhysicsState - '" + physicsName_ + "'";
}

template <class Scalar>
void PhysicsState<Scalar>::describe(
    Teuchos::FancyOStream& out,
    const Teuchos::EVerbosityLevel /* verbLevel */) const
{
  auto l_out = Teuchos::fancyOStream(out.getOStream());
  Teuchos::OSTab ostab(*l_out, 2, this->description());
  l_out->setOutputToRootOnly(0);

  *l_out << "\n--- " << this->description() << " ---" << std::endl;
}

}  // namespace Tempus
#endif  // Tempus_PhysicsState_impl_hpp
