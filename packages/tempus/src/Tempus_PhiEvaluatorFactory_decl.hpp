//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2026 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhiEvaluatorFactory_decl_hpp
#define Tempus_PhiEvaluatorFactory_decl_hpp

#include "Tempus_PhiEvaluator.hpp"

#include "Thyra_ModelEvaluator.hpp"

namespace Tempus {

/** \brief PhiEvaluator factory.
 *
 */
template <class Scalar>
class PhiEvaluatorFactory {
 public:
  /// Constructor
  PhiEvaluatorFactory() {}

  /// Destructor
  virtual ~PhiEvaluatorFactory() {}

  /// \name PhiEvaluator constructors
  //@{
  /// Create PhiEvaluator from PhiEvaluator type.
  Teuchos::RCP<Tempus::PhiEvaluator<Scalar> > createPhiEvaluator(
      std::string phiEvaluatorType = "PFD",
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model =
          Teuchos::null);

  /// Create PhiEvaluator from a ParameterList.
  Teuchos::RCP<Tempus::PhiEvaluator<Scalar> > createPhiEvaluator(
      Teuchos::RCP<Teuchos::ParameterList> phiEvaluatorPL,
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model =
          Teuchos::null);

 private:
  /// PhiEvaluator Factory.
  Teuchos::RCP<Tempus::PhiEvaluator<Scalar> > createPhiEvaluator(
      std::string phiEvaluatorType, Teuchos::RCP<Teuchos::ParameterList> phiEvaluatorPL,
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model);
};


}  // namespace Tempus

#endif  // Tempus_PhiEvaluatorFactory_decl_hpp
