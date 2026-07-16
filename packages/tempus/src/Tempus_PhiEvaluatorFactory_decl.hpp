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

/** \brief Creates supported PhiEvaluator implementations.
 *
 * The recognized type strings are "PFD", "Leja", and "Taylor".  The
 * type overload defaults an empty string to "PFD"; the parameter-list
 * overload reads "PhiEvaluator Type", also defaulting to "PFD".
 *
 * @tparam Scalar Scalar type of the evaluator and model.
 */
template <class Scalar>
class PhiEvaluatorFactory {
 public:
  /** \brief Constructs an empty factory. */
  PhiEvaluatorFactory() {}

  /** \brief Destroys the factory. */
  virtual ~PhiEvaluatorFactory() {}

  /// \name PhiEvaluator constructors
  //@{
  /** \brief Creates an evaluator selected by a type string.
   *
   * @param phiEvaluatorType std::string type selector: "PFD", "Leja", or
   * "Taylor".  An empty string selects "PFD".
   * @param model Const Teuchos::RCP to a Thyra::ModelEvaluator<Scalar>.
   * Currently ignored by this implementation.
   * @return Teuchos::RCP owning the selected PhiEvaluator<Scalar>.
   */
  Teuchos::RCP<Tempus::PhiEvaluator<Scalar> > createPhiEvaluator(
      std::string phiEvaluatorType = "PFD",
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model =
          Teuchos::null);

  /** \brief Creates an evaluator configured from a parameter list.
   *
   * The "PhiEvaluator Type" entry selects "PFD" when absent or when
   * @p phiEvaluatorPL is Teuchos::null.  The list is forwarded to the
   * selected evaluator when nonnull.
   *
   * @param phiEvaluatorPL Teuchos::RCP to the evaluator ParameterList, or
   * Teuchos::null.
   * @param model Const Teuchos::RCP to a Thyra::ModelEvaluator<Scalar>.
   * @return Teuchos::RCP owning the selected PhiEvaluator<Scalar>.
   */
  Teuchos::RCP<Tempus::PhiEvaluator<Scalar> > createPhiEvaluator(
      Teuchos::RCP<Teuchos::ParameterList> phiEvaluatorPL,
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model =
          Teuchos::null);

 private:
  /** \brief Dispatches to the selected nonmember evaluator constructor.
   *
   * @param phiEvaluatorType std::string type selector.
   * @param phiEvaluatorPL Teuchos::RCP to the optional configuration list.
   * @param model Const model RCP.
   * @return Teuchos::RCP owning the selected evaluator.
   */
  Teuchos::RCP<Tempus::PhiEvaluator<Scalar> > createPhiEvaluator(
      std::string phiEvaluatorType, Teuchos::RCP<Teuchos::ParameterList> phiEvaluatorPL,
      const Teuchos::RCP<const Thyra::ModelEvaluator<Scalar> >& model);
};


}  // namespace Tempus

#endif  // Tempus_PhiEvaluatorFactory_decl_hpp
