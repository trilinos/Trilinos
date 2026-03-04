//@HEADER
// *****************************************************************************
// TODO
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhiEvaluatorLeja_decl_hpp
#define Tempus_PhiEvaluatorLeja_decl_hpp

#include "Tempus_PhiEvaluator.hpp"

namespace Tempus {

/** \brief PhiEvaluatorLeja uses a Leja point method to evaluate the phi-function vector product
 *
 *  \f$[x = \varphi_k(J) b]\f$, where 
 *
 *   - b is a right hand side vector
 *   - J is a linear operator
 */
template <class Scalar>
class PhiEvaluatorLeja
  : virtual public PhiEvaluator<Scalar> {
 public:
  /// Inherit Contructor
  using PhiEvaluator<Scalar>::PhiEvaluator;

  /// \name Basic PhiEvaluatorLeja Methods
  //@{

  /// Return a valid ParameterList with current settings.
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const override;

 protected:
  Thyra::SolveStatus<Scalar> computeLinOpPhi(const int phi_order,
					     const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> L,
					     const Teuchos::RCP<Thyra::VectorBase<Scalar>> v) override;
};

/// Nonmember constructor from a ParameterList
template <class Scalar>
Teuchos::RCP<PhiEvaluatorLeja<Scalar> > createPhiEvaluatorLeja(
    Teuchos::RCP<Teuchos::ParameterList> pList);

}  // namespace Tempus

#endif  // Tempus_PhiEvaluatorLeja_decl_hpp
