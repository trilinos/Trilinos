//@HEADER
// *****************************************************************************
// TODO
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhiEvaluatorLeja_decl_hpp
#define Tempus_PhiEvaluatorLeja_decl_hpp

#include "Tempus_PhiEvaluator.hpp"

namespace Tempus {

enum LpType {
  LPREAL,
  LPCPLX,
  LPCONJ,
};

/*
 * Leja point container
 */
struct LejaPoint {
  std::complex<double> lp;
  LpType lpt;

  std::vector<std::complex<double>> get()
  {
    std::vector<std::complex<double>> out;
    switch (this->lpt) {
      case LPREAL:
        out.push_back(this->lp);
        break;
      case LPCONJ:
        out.push_back(this->lp);
        out.push_back(this->lp * std::complex(0.0, -1.0));
        break;
      case LPCPLX:
        out.push_back(this->lp);
        break;
    };
    return out;
  }
};

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
  PhiEvaluatorLeja<Scalar>(std::string name) : PhiEvaluator<Scalar>(name), maxLejaOrder_(0)
  { }
  PhiEvaluatorLeja<Scalar>() : PhiEvaluatorLeja<Scalar>("Phi Evaluator")
  { }

  /// \name Basic PhiEvaluatorLeja Methods
  //@{

  /// Set the linearization point for the Jacobian calculation
  // TODO: overwrite this, and compute optimal ellipse size in the future
  //void setLinearizationPoint(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &inArgs) overwrite;

  /// Return a valid ParameterList with current settings.
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const override;

  /// Set the parameters from a ParameterList
  void setPhiEvaluatorValues(Teuchos::RCP<Teuchos::ParameterList> pl) override;

  /// compute the shift and scale parameters from ellipse a,b,c bounds
  std::tuple<Scalar, Scalar> getShiftScale();

  /// Set the polynomial expansion order
  void setExpansionOrder(int order) { expansionOrder_ = order; }

  /// Get the polynomial expansion order
  int getExpansionOrder() const { return expansionOrder_; }

 protected:
  Thyra::SolveStatus<Scalar> computeLinOpPhi(const int phi_order,
					     const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> L,
					     const Teuchos::Ptr<Thyra::VectorBase<Scalar>> v) override;
 private:
  int maxLejaOrder_;
  int expansionOrder_;
  double leja_tol_;
  double leja_a_;
  double leja_b_;
  double leja_c_;

  Teuchos::Array<LejaPoint> lp_base_;
  Teuchos::Array<std::complex<double>> lp_dd_;

  /// TODO: compute the divided differences
  Teuchos::Array<std::complex<double>> updateDividedDiffs(const int phi_order);
};

/// Nonmember constructor from a ParameterList
template <class Scalar>
Teuchos::RCP<PhiEvaluatorLeja<Scalar> > createPhiEvaluatorLeja(
    Teuchos::RCP<Teuchos::ParameterList> pList);

}  // namespace Tempus

#endif  // Tempus_PhiEvaluatorLeja_decl_hpp
