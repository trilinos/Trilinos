//@HEADER
// *****************************************************************************
// TODO
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhiEvaluatorTaylor_decl_hpp
#define Tempus_PhiEvaluatorTaylor_decl_hpp

#include "Tempus_PhiEvaluator.hpp"
#include "Thyra_ProductVectorBase.hpp"


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
  using PhiEvaluator<Scalar>::PhiEvaluator;

  /// \name Basic PhiEvaluatorTaylor Methods

  /// Return a valid ParameterList with current settings.
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const override;

  /// Set the parameters from a ParameterList
  void setPhiEvaluatorValues(Teuchos::RCP<Teuchos::ParameterList> pl) override;

  void setTaylorExpansionOrder(int order) { taylorExpOrder_ = order; }
  int getTaylorExpansionOrder() const { return taylorExpOrder_; }

  void setLinearizationPoint(const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs) override;

  Thyra::SolveStatus<Scalar> computePhis(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> x,
					 Scalar cdt,
					 const std::vector<Teuchos::RCP<const Thyra::VectorBase<Scalar>>> rhs_B) override;

  Thyra::SolveStatus<Scalar> matrixExponential(const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> L,
					       const Teuchos::RCP<Thyra::VectorBase<Scalar>> v);

 private:
  mutable Teuchos::RCP<const Thyra::ModelEvaluatorBase::InArgs<Scalar>> inArgs_lin_;

  int taylorExpOrder_;

  Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> Atilde_;
  Teuchos::RCP<Thyra::VectorBase<Scalar>> v_;
  Teuchos::RCP<const Thyra::VectorBase<Scalar>> matExp_v_;
};

/// Nonmember constructor from a ParameterList
template <class Scalar>
Teuchos::RCP<PhiEvaluatorTaylor<Scalar>> createPhiEvaluatorTaylor(
    Teuchos::RCP<Teuchos::ParameterList> pList);

}  // namespace Tempus

#endif  // Tempus_PhiEvaluatorTaylor_decl_hpp
