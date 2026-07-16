//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2026 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhiEvaluatorPFD_decl_hpp
#define Tempus_PhiEvaluatorPFD_decl_hpp

#include "Tempus_PhiEvaluator.hpp"

namespace Tempus {

/** \brief Partial-fraction evaluator with a currently limited CN path.
 *
 * The implemented path evaluates only \f$\varphi_1\f$ by solving the
 * Crank-Nicolson form
 * \f$(M + (cdt/2)J_{\mathrm{impl}})x = M b\f$ using the prepared mass
 * matrix \f$M\f$ and implicit residual Jacobian \f$J_{\mathrm{impl}}\f$. It
 * does not implement the direct
 * linear-operator interface or other phi orders.  The "PFD Method" entry
 * accepts "CN"; other values currently issue a warning and use CN.
 *
 * @tparam Scalar Scalar type of the Thyra vectors and operators.
 */
template <class Scalar>
class PhiEvaluatorPFD
  : virtual public PhiEvaluator<Scalar> {
 public:
   /** \brief Inherits PhiEvaluator<Scalar> constructors. */
  using PhiEvaluator<Scalar>::PhiEvaluator;

  /// \name Basic PhiEvaluatorPFD Methods

   /** \brief Returns valid parameters, including the "PFD Method" default "CN". */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const override;

   /** \brief Applies base evaluator parameters and the PFD method selection.
    *
    * Only "CN" is implemented.  "Lump Mass Matrix" is disabled when set.
    *
    * @param pl Teuchos::RCP to the mutable evaluator ParameterList.
    */
  void setPhiEvaluatorValues(Teuchos::RCP<Teuchos::ParameterList> pl);

   /** \brief Computes the supported \f$\varphi_1\f$ action with CN.
    *
    * @pre phi_order is 1 and the evaluator has an initialized mass matrix and
    * Jacobian.
    * @param x Nonnull Teuchos::Ptr to the output Thyra::VectorBase<Scalar>.
    * @param phi_order int phi-function index; only 1 is accepted.
    * @param cdt Scalar time-scaled coefficient used as \f$cdt/2\f$ in the
    * Jacobian term.
    * @param Mrhs_b Const Teuchos::RCP to the mass-weighted right-hand-side
    * Thyra::VectorBase<Scalar>.
    * @return Thyra::SolveStatus<Scalar> from the mass-plus-Jacobian solve.
    */
  Thyra::SolveStatus<Scalar> computePhi(
    const Teuchos::Ptr<Thyra::VectorBase<Scalar>> x,
    const int phi_order, Scalar cdt,
    const Teuchos::RCP<const Thyra::VectorBase<Scalar>> &Mrhs_b) override;

   /** \brief Computes the supported \f$\varphi_1\f$ action from a RHS array.
    *
    * @pre Mrhs_B has length 2, Mrhs_B[0] is Teuchos::null, and Mrhs_B[1] is
    * the nonnull mass-weighted \f$\varphi_1\f$ right-hand side.
    * @param x Nonnull Teuchos::Ptr to the output Thyra::VectorBase<Scalar>.
    * @param cdt Scalar time-scaled coefficient used by the CN solve.
    * @param Mrhs_B Teuchos::ArrayView of exactly two vector RCPs, indexed by
    * phi order.
    * @return Thyra::SolveStatus<Scalar> from computePhi().
    */
  Thyra::SolveStatus<Scalar> computePhis(
    const Teuchos::Ptr<Thyra::VectorBase<Scalar>> x,
    const Scalar cdt,
    const Teuchos::ArrayView<const Teuchos::RCP<const Thyra::VectorBase<Scalar>>> &Mrhs_B) override;

  protected:
   /** \brief Rejects direct linear-operator evaluation.
    *
    * This override always throws std::invalid_argument.
    */
  Thyra::SolveStatus<Scalar> computeLinOpPhi(
    const int phi_order,
    const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> L,
    const Teuchos::Ptr<Thyra::VectorBase<Scalar>> v,
    const Scalar cdt=1.0
    ) override
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true,
                               std::invalid_argument,
                               "PhiEvaluatorPFD<Scalar>::computeLinOpPhi is not implemented.");
    return Thyra::SolveStatus<Scalar>();
  }
};

/** \brief Creates a PFD evaluator and applies an optional parameter list.
 *
 * @tparam Scalar Scalar type of the evaluator.
 * @param pList Teuchos::RCP to configuration, or Teuchos::null.
 * @return Teuchos::RCP owning the new PhiEvaluatorPFD<Scalar>.
 */
template <class Scalar>
Teuchos::RCP<PhiEvaluatorPFD<Scalar>> createPhiEvaluatorPFD(
    Teuchos::RCP<Teuchos::ParameterList> pList);

}  // namespace Tempus

#endif  // Tempus_PhiEvaluatorPFD_decl_hpp
