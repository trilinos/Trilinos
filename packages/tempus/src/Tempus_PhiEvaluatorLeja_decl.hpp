//@HEADER
// *****************************************************************************
// TODO
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhiEvaluatorLeja_decl_hpp
#define Tempus_PhiEvaluatorLeja_decl_hpp

#include "Tempus_PhiEvaluator.hpp"

#include <complex>

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
  // TODO: how does std::complex<double> interact with Scalar. Leja points are always complex.
  //       Template this on the Scalar tye and get the appropriate complex Scalar types?
  std::complex<double> lp;
  LpType lpt;

  std::vector<std::complex<double>> get()
  {
    std::vector<std::complex<double>> out;
    switch (this->lpt) {
      case LPREAL:
        out.push_back(std::complex(std::real(this->lp), 0.));
        break;
      case LPCONJ:
        out.push_back(this->lp);
        out.push_back(std::conj(this->lp));
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

  /// Update the Leja ellipse parameters
  void setLejaEllipse(Scalar a, Scalar b, Scalar c);

  /// Get shifted and scaled leja point (z_i)
  LejaPoint getLpSc(uint i);

  /// Compute divided differences
  Teuchos::ArrayRCP<std::complex<double>> getDividedDiffs(const int phi_order, const Scalar cdt, const int exp_order);

  /// Get the polynomial expansion order
  int getExpansionOrder() const { return expansionOrder_; }

 protected:
  Thyra::SolveStatus<Scalar> computeLinOpPhi(const int phi_order,
					     const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> L,
					     const Teuchos::Ptr<Thyra::VectorBase<Scalar>> v,
               const Scalar cdt=1.0
					     ) override;
 private:
  int maxLejaOrder_;
  int expansionOrder_;
  int ddMethod_;
  Scalar leja_tol_;
  Scalar leja_a_;
  Scalar leja_b_;
  Scalar leja_c_;

  /// Initialize the Leja points
  void initLejaPointsBase();

  /// Computes the divided differences via taylor series
  Teuchos::ArrayRCP<std::complex<double>> getDividedDiffsTS(const int phi_order, const Scalar cdt, const int exp_order);

  /// Computes the divided differences via recurrence relation
  Teuchos::ArrayRCP<std::complex<double>> getDividedDiffsRC(const int phi_order, const Scalar cdt, const int exp_order);

  /// Storage for the base Leja points
  Teuchos::ArrayRCP<LejaPoint> lejaPointsBase_;

  /// Storage for the Leja points
  Teuchos::ArrayRCP<LejaPoint> lp_;
};

/// Nonmember constructor from a ParameterList
template <class Scalar>
Teuchos::RCP<PhiEvaluatorLeja<Scalar> > createPhiEvaluatorLeja(
    Teuchos::RCP<Teuchos::ParameterList> pList);

}  // namespace Tempus

#endif  // Tempus_PhiEvaluatorLeja_decl_hpp
