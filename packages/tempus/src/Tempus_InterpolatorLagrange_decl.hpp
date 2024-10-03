//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_InterpolatorLagrange_decl_hpp
#define Tempus_InterpolatorLagrange_decl_hpp

#include "Tempus_config.hpp"
#include "Tempus_Interpolator.hpp"

namespace Tempus {

/** \brief Concrete implemenation of <tt>Interpolator</tt> that does
 * simple lagrange interpolation.
 */
template <class Scalar>
class InterpolatorLagrange : virtual public Interpolator<Scalar> {
 public:
  /// Contructor
  InterpolatorLagrange() : order_(0) {}

  /// Destructor
  ~InterpolatorLagrange() {}

  /// \name Overridden from Tempus::Interpolator
  //@{

  /// Store pointer to interpolation nodes
  void setNodes(
      const Teuchos::RCP<
          const std::vector<Teuchos::RCP<SolutionState<Scalar> > > >& nodes)
  {
    nodes_ = nodes;
  }

  /// Perform an interpolation
  void interpolate(const Scalar& t, SolutionState<Scalar>* state_out) const;

  /// Return the order of the interpolation
  int order() const { return order_; }

  //@}

  /// \name Overridden from Teuchos::Describable
  //@{
  std::string description() const { return "Tempus::InterpolatorLagrange"; }
  void describe(Teuchos::FancyOStream& out,
                const Teuchos::EVerbosityLevel /* verbLevel */) const
  {
    out.setOutputToRootOnly(0);
    out << description() << "::describe" << std::endl;
  }
  //@}

  /// \name Overridden from Teuchos::ParameterListAcceptor
  //@{
  void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList);
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
  //@}

 private:
  void lagrange(const int p, const Scalar& t,
                SolutionState<Scalar>* state_out) const;

  Teuchos::RCP<const std::vector<Teuchos::RCP<SolutionState<Scalar> > > >
      nodes_;
  Teuchos::RCP<Teuchos::ParameterList> pl_;
  int order_;
};

// Nonmember constructor
template <class Scalar>
Teuchos::RCP<InterpolatorLagrange<Scalar> > lagrangeInterpolator()
{
  return Teuchos::rcp(new InterpolatorLagrange<Scalar>);
}

}  // namespace Tempus

#endif  // Tempus_InterpolatorLagrange_decl_hpp
