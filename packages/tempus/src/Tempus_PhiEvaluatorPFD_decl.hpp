//@HEADER
// *****************************************************************************
// TODO
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhiEvaluatorPFD_decl_hpp
#define Tempus_PhiEvaluatorPFD_decl_hpp

#include "Tempus_PhiEvaluator.hpp"

namespace Tempus {

/** \brief PhiEvaluatorPFD uses a partial fraction decomposition to evaluate
 *
 *  \f$[x = \varphi_k(J) b]\f$, where
 *
 *   - b is a right hand side vector
 *   - J is a linear operator
 */
template <class Scalar>
class PhiEvaluatorPFD
  : virtual public PhiEvaluator<Scalar> {
 public:
  /// Inherit Contructor
  using PhiEvaluator<Scalar>::PhiEvaluator;

  /// \name Basic PhiEvaluatorPFD Methods

  /// Return a valid ParameterList with current settings.
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const override;

  /// Set the parameters from a ParameterList
  void setPhiEvaluatorValues(Teuchos::RCP<Teuchos::ParameterList> pl);

  /// compute the Phi_k function of cdt times Jacobian for right hand side rhs_b
  Thyra::SolveStatus<Scalar> computePhi(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> x,
					const int phi_order, Scalar cdt,
					const Teuchos::RCP<const Thyra::VectorBase<Scalar>> &Mrhs_b) override;

  /// compute the Phi_k function of cdt times Jacobian for a linear combination with right hand side vector rhs_B
  Thyra::SolveStatus<Scalar> computePhis(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> x,
                                         const Scalar cdt,
                                         const Teuchos::ArrayView<const Teuchos::RCP<const Thyra::VectorBase<Scalar>>> &Mrhs_B) override;

 protected:
  Thyra::SolveStatus<Scalar> computeLinOpPhi(const int phi_order,
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

/// Nonmember constructor from a ParameterList
template <class Scalar>
Teuchos::RCP<PhiEvaluatorPFD<Scalar>> createPhiEvaluatorPFD(
    Teuchos::RCP<Teuchos::ParameterList> pList);

}  // namespace Tempus

#endif  // Tempus_PhiEvaluatorPFD_decl_hpp
