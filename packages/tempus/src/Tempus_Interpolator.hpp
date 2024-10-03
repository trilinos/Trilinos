//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_Interpolator_hpp
#define Tempus_Interpolator_hpp

// Teuchos
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_ParameterListAcceptorDefaultBase.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"

// Tempus
#include "Tempus_config.hpp"
#include "Tempus_SolutionState.hpp"

#include <vector>

namespace Tempus {

/** \brief Base strategy class for interpolation functionality.
 */
template <class Scalar>
class Interpolator
  : virtual public Teuchos::Describable,
    virtual public Teuchos::ParameterListAcceptor,
    virtual public Teuchos::VerboseObject<Interpolator<Scalar> > {
 public:
  /** \brief Store pointer to interpolation nodes
   *
   * This function represent a persisting relationship between the
   * interpolation nodes and the interpolator.  For simple interpolators like
   * linear and Hermite, this is not needed, but for interpolators like cubic
   * splines where there is some computational work in assembling the
   * interpolant, this is important.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>nodes</tt> must have unique time values and be sorted in ascending
   * time order
   * </ul>
   *
   * <b>Postconditions:</b><ul>
   * <li>If this function is called twice and <tt>nodes</tt> is a different
   * pointer than was previously called, then it is possible that the
   * interpolant will be recomputed when <tt>interpolate</tt> is next called.
   * </ul>
   */
  virtual void setNodes(
      const Teuchos::RCP<
          const std::vector<Teuchos::RCP<SolutionState<Scalar> > > >&
          nodes) = 0;

  /** \brief Perform an interpolation.
   */
  virtual void interpolate(const Scalar& t,
                           SolutionState<Scalar>* state_out) const = 0;

  /** \brief Return the order of the interpolation.
   */
  virtual int order() const = 0;
};

/// Nonmember functions
template <class Scalar>
void interpolate(const Interpolator<Scalar>& interpolator, const Scalar& t,
                 SolutionState<Scalar>* state_out)
{
  interpolator.interpolate(t, state_out);
}

/// Nonmember functions
template <class Scalar>
void interpolate(
    Interpolator<Scalar>& interpolator,
    const Teuchos::RCP<
        const std::vector<Teuchos::RCP<SolutionState<Scalar> > > >& nodes,
    const Scalar& t, SolutionState<Scalar>* state_out)
{
  interpolator.setNodes(nodes);
  interpolator.interpolate(t, state_out);
}

/// Helper function for comparing times
template <typename Scalar>
bool floating_compare_equals(const Scalar& a, const Scalar& b,
                             const Scalar& scale)
{
  typedef Teuchos::ScalarTraits<Scalar> ST;
  typedef typename ST::magnitudeType mag_type;

  const mag_type tol  = ST::magnitude(scale) * ST::eps() * mag_type(1000.0);
  const mag_type diff = ST::magnitude(a - b);

  return diff <= tol || diff <= ST::sfmin();
}

}  // namespace Tempus

#endif  // Tempus_Interpolator_hpp
