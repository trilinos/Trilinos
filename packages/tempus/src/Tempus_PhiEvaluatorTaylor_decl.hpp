//@HEADER
// *****************************************************************************
// TODO
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhiEvaluatorTaylor_decl_hpp
#define Tempus_PhiEvaluatorTaylor_decl_hpp

#include "Tempus_PhiEvaluator.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Teuchos_TimeMonitor.hpp"


namespace Tempus {

/** \brief PhiEvaluatorTaylor uses a Taylor expansion to compute
 *
 *  \f$[x = \varphi_k(J) b]\f$, where
 *
 *   - b is a right hand side vector
 *   - J is a linear operator
 */
template <class Scalar>
class PhiEvaluatorTaylor
  : virtual public PhiEvaluator<Scalar> {
 public:
  /// Inherit Contructor
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

  /// Return a valid ParameterList with current settings.
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const override;

  /// Set the parameters from a ParameterList
  void setPhiEvaluatorValues(Teuchos::RCP<Teuchos::ParameterList> pl) override;

  /// Set the polynomial expansion order
  void setExpansionOrder(const int expansionOrder) { expansionOrder_ = expansionOrder; }

  /// Get the polynomial expansion order
  int getExpansionOrder() const { return expansionOrder_; }

 protected:
  Thyra::SolveStatus<Scalar> computeLinOpPhi(const int phi_order,
               const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> L,
               const Teuchos::Ptr<Thyra::VectorBase<Scalar>> v,
               const Scalar cdt=1.0
               ) override;

 private:
  int expansionOrder_;

#ifdef TEMPUS_TEUCHOS_TIME_MONITOR
  Teuchos::RCP<Teuchos::Time> timerLinOp_, timerPhi_;
#endif
};

/// Nonmember constructor from a ParameterList
template <class Scalar>
Teuchos::RCP<PhiEvaluatorTaylor<Scalar>> createPhiEvaluatorTaylor(
    Teuchos::RCP<Teuchos::ParameterList> pList);

}  // namespace Tempus

#endif  // Tempus_PhiEvaluatorTaylor_decl_hpp
