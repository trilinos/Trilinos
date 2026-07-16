//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2026 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhiEvaluatorTaylor_decl_hpp
#define Tempus_PhiEvaluatorTaylor_decl_hpp

#include "Tempus_PhiEvaluator.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Teuchos_TimeMonitor.hpp"


namespace Tempus {

/** \brief Evaluates phi functions with a truncated Taylor series.
 *
 * For a linear operator \f$L\f$, this evaluator forms
 * \f$\varphi_p(L)v = \sum_{j=0}^{N} L^j v/(j+p)!\f$, where \f$p\f$ is the
 * phi order and \f$N\f$ is "Expansion Order".  Single-RHS phi evaluations
 * use this series directly.  Multi-RHS evaluations use the base-class
 * augmented-operator extension and therefore evaluate \f$\varphi_0\f$ of
 * that extension.
 *
 * @tparam Scalar Scalar type of the Thyra vectors and operators.
 */
template <class Scalar>
class PhiEvaluatorTaylor
  : virtual public PhiEvaluator<Scalar> {
 public:
   /** \brief Constructs an evaluator with the supplied descriptive name.
    * @param name std::string stored as the evaluator name.
    */
  PhiEvaluatorTaylor<Scalar>(std::string name)
    : PhiEvaluator<Scalar>(name),
      expansionOrder_(0)
  {
#ifdef TEMPUS_TEUCHOS_TIME_MONITOR
    std::stringstream ss;
    ss << "Tempus::" << name ;

    std::string phiLabel = ss.str() + ": PhiEval";
    timerPhi_ = Teuchos::TimeMonitor::getNewCounter(phiLabel);

    std::string linOpLabel = ss.str() + ": LinOp";
    timerLinOp_ = Teuchos::TimeMonitor::getNewCounter(linOpLabel);
#endif
  }
  PhiEvaluatorTaylor<Scalar>() : PhiEvaluatorTaylor<Scalar>("PhiEvaluatorTaylor")
  { }

   /// \name Basic PhiEvaluatorTaylor Methods

  std::string description() const;
  void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const;

   /** \brief Returns valid parameters including defaults. */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const override;

   /** \brief Applies base parameters and "Expansion Order".
    *
    * @param pl Teuchos::RCP to the mutable evaluator ParameterList.
    */
  void setPhiEvaluatorValues(Teuchos::RCP<Teuchos::ParameterList> pl) override;

   /** \brief Sets the int Taylor degree \f$N > 0\f$ used for new evaluations.
    * @param expansionOrder Number of operator-power updates after the initial
    * \f$j=0\f$ term.
    */
  void setExpansionOrder(const int expansionOrder) { expansionOrder_ = expansionOrder; }

   /** \brief Returns the int Taylor degree \f$N>0\f$. */
  int getExpansionOrder() const { return expansionOrder_; }

  protected:
   /** \brief Applies the truncated Taylor series in place.
    *
    * @param phi_order Nonnegative int phi-function index \f$p\f$.
    * @param L Const Teuchos::RCP to the already time-scaled
    * Thyra::LinearOpBase<Scalar>.
    * @param v Nonnull Teuchos::Ptr to the input vector, overwritten by the
    * series result.
    * @param cdt Scalar retained by the virtual interface; this implementation
    * uses the scaling already contained in @p L.
    * @return Thyra::SolveStatus<Scalar> whose achieved tolerance is the final
    * update infinity norm.
    */
  Thyra::SolveStatus<Scalar> computeLinOpPhi(const int phi_order,
               const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> L,
               const Teuchos::Ptr<Thyra::VectorBase<Scalar>> v,
               const Scalar cdt=1.0
               ) override;

  private:
  /** \brief int Taylor degree \f$N>0\f$ configured by "Expansion Order". */
  int expansionOrder_;

#ifdef TEMPUS_TEUCHOS_TIME_MONITOR
  Teuchos::RCP<Teuchos::Time> timerLinOp_, timerPhi_;
#endif
};

/** \brief Creates a Taylor evaluator and applies an optional parameter list.
 *
 * @tparam Scalar Scalar type of the evaluator.
 * @param pList Teuchos::RCP to configuration, or Teuchos::null.
 * @return Teuchos::RCP owning the new PhiEvaluatorTaylor<Scalar>.
 */
template <class Scalar>
Teuchos::RCP<PhiEvaluatorTaylor<Scalar>> createPhiEvaluatorTaylor(
    Teuchos::RCP<Teuchos::ParameterList> pList);

}  // namespace Tempus

#endif  // Tempus_PhiEvaluatorTaylor_decl_hpp
