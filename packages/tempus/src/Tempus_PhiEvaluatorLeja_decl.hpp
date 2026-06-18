//@HEADER
// *****************************************************************************
// TODO
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhiEvaluatorLeja_decl_hpp
#define Tempus_PhiEvaluatorLeja_decl_hpp

#include "Tempus_PhiEvaluator.hpp"
#include "Teuchos_TimeMonitor.hpp"

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
  // A Leja point is based on a point on the complex plane.
  // It is always a complex double, regardless of Scalar type.
  // Supporting complex float might be possible, but may not be beneficial, since operations on Leja points are cheap.
  std::complex<double> lp;

  // LPREAL means that lp is supposed to be real, any imaginary component must be disregarded.
  // LPCPLX means that lp is complex.
  // LPCONJ means that lp is complex, and represents a conjugate pair, both lp and std::conj(lp).
  LpType lpt;

  // get the list of sanitized and expanded Leja points
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
  PhiEvaluatorLeja<Scalar>(std::string name)
    : PhiEvaluator<Scalar>(name),
      ddMethod_(0),
      leja_sf_(1.0),
      leja_tol_(1.0e-6),
      leja_a_(-1.0),
      leja_b_(0.0),
      leja_c_(1.0)
  {
#ifdef TEMPUS_TEUCHOS_TIME_MONITOR
    std::stringstream ss;
    ss << "Tempus::" << name;

    // set up timers for overall Phi-evaluation, divided differences, and Linop-vector applications
    std::string phiLabel = ss.str() + ": PhiEval";
    timerPhi_ = Teuchos::TimeMonitor::getNewCounter(phiLabel);
    std::string linOpLabel = ss.str() + ": LinOp";
    timerLinOp_ = Teuchos::TimeMonitor::getNewCounter(linOpLabel);
    std::string ddLabel = ss.str() + ": dd_phi";
    timerDD_ = Teuchos::TimeMonitor::getNewCounter(ddLabel);
#endif
  }
  PhiEvaluatorLeja<Scalar>() : PhiEvaluatorLeja<Scalar>("PhiEvaluatorLeja")
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

  /// Compute the scale, normalized shift and normalized anisotropy
  /// parameters from ellipse a,b,c bounds.
  constexpr std::tuple<double, double, double> getScaleFromBase();

  // transform base LejaPoint to shifted anisotropic transformed Leja point (not scaled)
  constexpr LejaPoint transformLejaPoint(const LejaPoint& lp_base, const std::tuple<const double, const double, const double>& scale_params);

  /// Set the polynomial expansion order
  void setExpansionOrder(const int expansionOrder);

  /// Update the Leja ellipse parameters a,b,c
  void setLejaEllipse(const double a, const double b, const double c);

  /// Get ellipse parameters a,b,c
  Teuchos::Tuple<double, 3> getLejaEllipse();

  /// Adaptively update the spectrum parameters a,b,c using Krylov-Schur
  void adaptEvaluator() override;

  /// Set the divided difference method
  void setDivideDifferenceMethod(const int ddMethod);

  /// Get shifted and scaled leja point (z_i)
  LejaPoint getLpSc(const int lp_idx);

  /// Compute divided differences
  Teuchos::ArrayRCP<Scalar> getDividedDiffs(const int phi_order, const Scalar cdt, const int lejaOrder);

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
  int ddMethod_;

  // those are read as double from the input file regardless of scalar type
  double leja_sf_;
  double leja_tol_;
  double leja_a_;
  double leja_b_;
  double leja_c_;

#ifdef TEMPUS_TEUCHOS_TIME_MONITOR
  Teuchos::RCP<Teuchos::Time> timerDD_, timerLinOp_, timerPhi_;
#endif

  /// allows to check at runtime for compelex Scalar type
  static const bool isComplex = Teuchos::ScalarTraits<Scalar>::isComplex;

  /// Initialize lejaOrder Leja points (for negative lejaOrder, deduce the number of points from expansionOrder_)
  void initLejaPointsBase(const int lejaOrder = -1);

  /// Computes the divided differences via taylor series
  Teuchos::ArrayRCP<Scalar> getDividedDiffsTS(const int phi_order, const Scalar cdt, const int lejaOrder);

  /// Computes the divided differences via taylor series
  Teuchos::ArrayRCP<Scalar> getDividedDiffsTSR(const int phi_order, const Scalar cdt, const int lejaOrder);

  /// Computes the divided differences via recurrence relation
  Teuchos::ArrayRCP<Scalar> getDividedDiffsRC(const int phi_order, const Scalar cdt, const int lejaOrder);

  /// Computes divided differences of phi_k via the Zivcovich (2019) dd_phi method.
  /// Ref: F. Zivcovich. "Fast and accurate computation of divided differences for analytic
  /// functions, with an application to the exponential function."
  /// Dolomites Research Notes on Approximation. 12. 2019.
  Teuchos::ArrayRCP<Scalar> getDividedDiffsPhi(const int phi_order, const Scalar cdt, const int lejaOrder);

  /// Storage for the base Leja points
  Teuchos::ArrayRCP<LejaPoint> lejaPointsBase_;
};

/// Nonmember constructor from a ParameterList
template <class Scalar>
Teuchos::RCP<PhiEvaluatorLeja<Scalar> > createPhiEvaluatorLeja(
    Teuchos::RCP<Teuchos::ParameterList> pList);

}  // namespace Tempus

#endif  // Tempus_PhiEvaluatorLeja_decl_hpp
