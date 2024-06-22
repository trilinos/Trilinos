//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhysicsState_hpp
#define Tempus_PhysicsState_hpp

#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include <string>
#include "Tempus_config.hpp"

namespace Tempus {

template <class Scalar>
/** \brief PhysicsState is a simple class to hold information about the physics.
 *
 * <b>Design Considerations</b>
 *   - The purpose of PhysicsState is to provide a means to store any physics
 *     information that is required to restart from a checkpoint.
 *   - The PhysicsState should be held by the SolutionState so that it has
 *     all the information needed to restart the solution.
 *   - Many Physics will simply need this base class, because
 *     they do not have any other additional state information.
 *   - PhysicsState can be inherited to expand the state information.
 *   - Examples of other information that could be included in derived
 *     PhysicsStates:
 *     - Specialized data structures (e.g., particle lists)
 *   - The base class currently has a name so the Physics can be checked to
 *     see if PhysicsState matches the current physics.
 */
class PhysicsState
  : public Teuchos::Describable,
    public Teuchos::VerboseObject<Tempus::PhysicsState<Scalar> > {
 public:
  /// Constructor
  PhysicsState(std::string pN = "Tempus::PhysicsState");

  /// Destructor
  virtual ~PhysicsState() {}

  /// Clone constructor
  virtual Teuchos::RCP<PhysicsState<Scalar> > clone() const;

  /// This is a deep copy
  virtual void copy(const Teuchos::RCP<const PhysicsState<Scalar> >& pS);

  /// Return name of PhysicsState
  virtual std::string getName() const;

  /// Set name of PhysicsState
  virtual void setName(std::string pN);

  /// \name Overridden from Teuchos::Describable
  //@{
  virtual std::string description() const;

  virtual void describe(Teuchos::FancyOStream& out,
                        const Teuchos::EVerbosityLevel verbLevel) const;
  //@}

 protected:
  std::string physicsName_;  ///< Name of the creating Physics.
};
}  // namespace Tempus
#endif  // Tempus_PhysicsState_hpp
