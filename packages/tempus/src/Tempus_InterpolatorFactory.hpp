//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_InterpolatorFactory_hpp
#define Tempus_InterpolatorFactory_hpp

#include "Teuchos_ParameterList.hpp"
#include "Tempus_config.hpp"
#include "Tempus_InterpolatorLagrange.hpp"

namespace Tempus {

/** \brief Interpolator factory.
 *
 * <b>Adding Interpoloators</b>
 *    -#
 */
template <class Scalar>
class InterpolatorFactory {
 public:
  /// Create default interpolator from interpolator type (e.g., "Linear").
  static Teuchos::RCP<Interpolator<Scalar> > createInterpolator(
      std::string interpolatorType = "")
  {
    if (interpolatorType == "") interpolatorType = "Lagrange";
    return createInterpolator(interpolatorType, Teuchos::null);
  }

  /// Create interpolator from ParameterList with its details.
  static Teuchos::RCP<Interpolator<Scalar> > createInterpolator(
      const Teuchos::RCP<Teuchos::ParameterList>& interpolatorPL)
  {
    std::string interpolatorType =
        interpolatorPL->get<std::string>("Interpolator Type", "Lagrange");
    return createInterpolator(interpolatorType, interpolatorPL);
  }

 private:
  /// Very simple factory method
  static Teuchos::RCP<Interpolator<Scalar> > createInterpolator(
      const std::string& interpolatorType,
      const Teuchos::RCP<Teuchos::ParameterList>& interpolatorPL)
  {
    using Teuchos::rcp;

    Teuchos::RCP<Interpolator<Scalar> > interpolator;
    if (interpolatorType == "Lagrange")
      interpolator = rcp(new InterpolatorLagrange<Scalar>);
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(
          true, std::logic_error,
          "Unknown 'Interpolator Type' = " << interpolatorType);
    }
    interpolator->setParameterList(interpolatorPL);

    return interpolator;
  }
};

}  // namespace Tempus
#endif  // Tempus_InterpolatorFactory_hpp
