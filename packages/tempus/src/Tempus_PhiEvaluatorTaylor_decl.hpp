//@HEADER
// *****************************************************************************
// TODO
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhiEvaluatorTaylor_decl_hpp
#define Tempus_PhiEvaluatorTaylor_decl_hpp

#include "Tempus_PhiEvaluator.hpp"

namespace Tempus {

/** \brief PhiEvaluatorTaylor uses a partial fraction decomposition to evaluate
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

  void setLinearizationPoint(const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs) override;

  Thyra::SolveStatus<Scalar> computePhi(const Teuchos::Ptr<Thyra::VectorBase<Scalar>> vphi,
					int k, Scalar cdt, const Teuchos::RCP<const Thyra::VectorBase<Scalar>> rhs_b) override;

 private:
  mutable Teuchos::RCP<const Thyra::ModelEvaluatorBase::InArgs<Scalar>> inArgs_lin_;
};

/// Nonmember constructor from a ParameterList
template <class Scalar>
Teuchos::RCP<PhiEvaluatorTaylor<Scalar>> createPhiEvaluatorTaylor(
    Teuchos::RCP<Teuchos::ParameterList> pList);

}  // namespace Tempus

#endif  // Tempus_PhiEvaluatorTaylor_decl_hpp
