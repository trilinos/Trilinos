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

/** \brief Partial-fraction evaluator.
 *
 * The PFD method evaluates \f$\varphi_s\f$ by solving a partial fraction form
 * \f$ r(J)f = p(J) f + \sum_{m=1}^{M'} \gamma_m v_m \f$,
 * where \f$ p \f$ is a small polynomial (usually zero or a constant) and
 * the updates (partial fractions) are given by
 * \f$ (\lambda_m M + J_{\mathrm{impl}}) v_m = Mf \f$
 * for \f$ m = 1,\ldots,M' \f$.
 *
 * @tparam Scalar Scalar type of the Thyra vectors and operators.
 */
template <class Scalar>
class PhiEvaluatorPFD
  : virtual public PhiEvaluator<Scalar> {
  /** \brief Structure to hold the coefficients of the partial fraction decomposition
   */
  struct RealPF {
    /// (real) location of the pole of the partial fraction
    double pole;
    /// coefficient of the partial fraction
    double coeff;
  };
  /// list of partial fractions of the expansion
  using PFList = std::vector<RealPF>;
  /// polynomial part (monomial coefficients, constant, linear, ....)
  using PFPolynomial = std::vector<double>;

  /// struct that defines a partial fraction decomposition
  struct PFDMethod {
    PFList pfs;
    PFPolynomial poly;
  };

 public:
  /** \brief Constructs an evaluator with the supplied descriptive name.
   * @param name std::string stored as the evaluator name.
   */
  PhiEvaluatorPFD(std::string name);
  PhiEvaluatorPFD() : PhiEvaluatorPFD("PhiEvaluatorPFD")
  {
  }

  /// \name Basic PhiEvaluatorPFD Methods

  /** \brief Returns valid parameters, including the "PFD Method" default "CN". */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const override;

  /** \brief Applies base evaluator parameters and the PFD method selection.
   *
   * @param pl Teuchos::RCP to the mutable evaluator ParameterList.
   */
  void setPhiEvaluatorValues(Teuchos::RCP<Teuchos::ParameterList> pl);

  /** \brief Computes the supported \f$\varphi_s\f$ action.
   *
   * @pre the evaluator has an initialized mass matrix and Jacobian.
   * @param phiv Nonnull Teuchos::Ptr to the output Thyra::VectorBase<Scalar>.
   * @param phi_order int phi-function index.
   * @param cdt Scalar time-scaled coefficient used in the Jacobian term.
   * @param Mrhs_b Const Teuchos::RCP to the mass-weighted right-hand-side
   * Thyra::VectorBase<Scalar>.
   * @return Thyra::SolveStatus<Scalar> from the mass-plus-Jacobian solve.
   */
  Thyra::SolveStatus<Scalar> computePhi(
      const Teuchos::Ptr<Thyra::VectorBase<Scalar>> phiv,
      const int phi_order, Scalar cdt,
      const Teuchos::RCP<const Thyra::VectorBase<Scalar>> &Mrhs_b) override;

  /** \brief Computes the \f$\varphi_k\f$ action from a RHS array.
   *
   * @pre Mrhs_B has length K, and Mrhs_B[k] is
   *      the mass-weighted \f$\varphi_k\f$ right-hand side or Teuchos::null.
   * @param phiv Nonnull Teuchos::Ptr to the output Thyra::VectorBase<Scalar>.
   * @param cdt Scalar time-scaled coefficient used by the CN solve.
   * @param Mrhs_B Teuchos::ArrayView of vector RCPs, indexed by phi order k.
   * @return Thyra::SolveStatus<Scalar>.
   */
  Thyra::SolveStatus<Scalar> computePhis(
      const Teuchos::Ptr<Thyra::VectorBase<Scalar>> phiv,
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
      const Scalar cdt = 1.0) override
  {
    TEUCHOS_TEST_FOR_EXCEPTION(true,
                               std::invalid_argument,
                               "PhiEvaluatorPFD<Scalar>::computeLinOpPhi is not implemented.");
    return Thyra::SolveStatus<Scalar>();
  }

 private:
  /// string of partial fraction expansion method
  std::string pfdMethod_ = "CN";

  /// maps from method string to partial fraction expansion
  std::map<std::string, PFDMethod> pfdMethods_;

#ifdef TEMPUS_TEUCHOS_TIME_MONITOR
  Teuchos::RCP<Teuchos::Time> timerPhi_;
  Teuchos::RCP<Teuchos::Time> timerLinOp_;
  Teuchos::RCP<Teuchos::Time> timerSolve_;
#endif
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
