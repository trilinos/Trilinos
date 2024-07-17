//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhysicsStateCounter_hpp
#define Tempus_PhysicsStateCounter_hpp

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include <string>
#include "Tempus_PhysicsState.hpp"

namespace Tempus_Test {

template <class Scalar>
/** \brief PhysicsStateCounter is a simple PhysicsState that counts steps.
 *
 */
class PhysicsStateCounter : virtual public Tempus::PhysicsState<Scalar> {
 public:
  /// Constructor
  PhysicsStateCounter(std::string pN = "Tempus::PhysicsStateCounter",
                      int pI         = 0)
    : Tempus::PhysicsState<Scalar>(pN), physicsCounter_(pI)
  {
  }

  /// Destructor
  virtual ~PhysicsStateCounter() {}

  /// Clone constructor
  virtual Teuchos::RCP<Tempus::PhysicsState<Scalar> > clone() const
  {
    Teuchos::RCP<PhysicsStateCounter<Scalar> > pSC =
        Teuchos::rcp(new PhysicsStateCounter<Scalar>(this->physicsName_,
                                                     this->physicsCounter_));
    return pSC;
  }

  // using Tempus::PhysicsState<Scalar>::copy;
  /// This is a deep copy
  virtual void copy(const Teuchos::RCP<const Tempus::PhysicsState<Scalar> >& pS)
  {
    Teuchos::RCP<const PhysicsStateCounter<Scalar> > pSC =
        Teuchos::rcp_dynamic_cast<const PhysicsStateCounter<Scalar> >(pS);

    this->physicsName_    = pSC->getName();
    this->physicsCounter_ = pSC->getCounter();
  }

  /// Return counter of PhysicsStateCounter
  virtual int getCounter() const { return physicsCounter_; }

  /// Set counter of PhysicsStateCounter
  virtual void setCounter(int counter) { physicsCounter_ = counter; }

  /// \name Overridden from Teuchos::Describable
  //@{
  virtual void describe(Teuchos::FancyOStream& out,
                        const Teuchos::EVerbosityLevel verbLevel) const
  {
    out << this->description() << "::describe" << std::endl
        << "  physicsName      = " << this->physicsName_ << std::endl
        << "  physicsCounter = " << physicsCounter_ << std::endl;
  }
  //@}

 protected:
  int physicsCounter_;  ///< Counter for steps
};
}  // namespace Tempus_Test
#endif  // Tempus_PhysicsStateCounter_hpp
