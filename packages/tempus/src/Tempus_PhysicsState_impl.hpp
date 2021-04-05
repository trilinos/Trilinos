// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Tempus_PhysicsState_impl_hpp
#define Tempus_PhysicsState_impl_hpp

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include <string>


namespace Tempus {

template<class Scalar>
PhysicsState<Scalar>::PhysicsState(std::string pN)
  : physicsName_(pN)
{}

template<class Scalar>
Teuchos::RCP<PhysicsState<Scalar> >
PhysicsState<Scalar>::clone() const
{
  Teuchos::RCP<PhysicsState<Scalar> > pSB =
    Teuchos::rcp(new PhysicsState<Scalar> ());

  pSB->setName(this->getName());

  return pSB;
}

template<class Scalar>
void PhysicsState<Scalar>::copy(
  const Teuchos::RCP<const PhysicsState<Scalar> >& pS)
{
  physicsName_ = pS->getName();
}

template<class Scalar>
std::string PhysicsState<Scalar>::getName() const
{
  return physicsName_;
}

template<class Scalar>
void PhysicsState<Scalar>::setName(std::string pN)
{
  physicsName_ = pN;
}

template<class Scalar>
std::string PhysicsState<Scalar>::description() const
{
  return physicsName_;
}

template<class Scalar>
void PhysicsState<Scalar>::describe(
  Teuchos::FancyOStream        & in_out,
  const Teuchos::EVerbosityLevel /* verbLevel */) const
{
  auto out = Teuchos::fancyOStream( in_out.getOStream() );
  out->setOutputToRootOnly(0);
  *out << description() << "::describe" << std::endl
      << "  physicsName   = " << physicsName_ << std::endl;
}


} // namespace Tempus
#endif // Tempus_PhysicsState_impl_hpp
