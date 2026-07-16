//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2026 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_PhiEvaluatorLeja_decl_hpp
#define Tempus_PhiEvaluatorLeja_decl_hpp

#include "Tempus_PhiEvaluator.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include <complex>

namespace Tempus {

/** \brief Interpretation of a stored LejaPoint complex value. */
enum LpType {
  /** \brief Real point; discard any stored imaginary component. */
  LPREAL,
  /** \brief One non-conjugate complex point. Currently not supported by the implementation. */
  LPCPLX,
  /** \brief Upper-half-plane point representing it and its conjugate. */
  LPCONJ,
};

/** \brief Stores one base or transformed Leja interpolation point.
 *
 * The point is always std::complex<double>, independent of Scalar.
 * LPCONJ stores one upper-half-plane point but expands to two interpolation
 * points; LPREAL and LPCPLX expand to one point.
 */
struct LejaPoint {
  /** \brief std::complex<double> coordinate of this interpolation point.
   *
   * A Leja point is based on a point on the complex plane.
   * It is always a complex double, regardless of Scalar type.
   * Supporting complex float might be possible, but may not be beneficial, since operations on Leja points are cheap.
   */
  std::complex<double> lp;

  /** \brief LpType controlling real sanitization or conjugate expansion.
   *
   * LPREAL means that lp is supposed to be real, any imaginary component must be disregarded.
   * LPCPLX means that lp is complex.
   * LPCONJ means that lp is complex, and represents a conjugate pair, both lp and std::conj(lp).
   */
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

/** \brief Evaluates matrix exponentials using a Leja interpolation polynomial.
 *
 * Base Leja points are shifted and anisotropically transformed using ellipse
 * parameters whose real-axis extrema are \f$a\f$ and \f$b\f$ and whose
 * imaginary semi-axis is \f$c\f$.  These parameters are selected from
 * estimated spectral bounds or specified manually.  The direct
 * operator path currently supports only
 * \f$\varphi_0(L)v = \exp(L)v\f$.  Phi functions of positive order are
 * handled through the base-class augmented-operator extension.
 *
 * "Leja DD Method" chooses divided differences: 0 recurrence, 1 complex
 * Taylor scaling-and-squaring, 2 dd_phi H-factorization, or 3 real Taylor
 * scaling-and-squaring.  "Leja Ellipse Safety Factor" affects adaptation by
 * multiplying the estimated \f$a\f$ and \f$c\f$ values; \f$b\f$ is unchanged.
 *
 * @tparam Scalar Scalar type of the Thyra vectors and operators.
 */
template <class Scalar>
class PhiEvaluatorLeja
  : virtual public PhiEvaluator<Scalar> {
 public:
   /** \brief Constructs an evaluator with the supplied descriptive name.
    * @param name std::string stored as the evaluator name.
    */
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

  std::string description() const;
  void describe(Teuchos::FancyOStream& out, const Teuchos::EVerbosityLevel verbLevel) const;

   /** \brief Returns valid Leja parameters and their defaults. */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const override;

   /** \brief Applies Leja configuration values from a parameter list.
    *
    * "Expansion Order" sets polynomial degree, "leja_tol" terminates on the
    * update infinity norm, and "leja_a", "leja_b", and "leja_c" define the
    * initial ellipse.  "Leja DD Method" selects the divided-difference
    * algorithm.
    *
    * @param pl Teuchos::RCP to the mutable evaluator ParameterList.
    */
  void setPhiEvaluatorValues(Teuchos::RCP<Teuchos::ParameterList> pl) override;

   /** \brief Computes transform parameters from the current ellipse.
    * @return std::tuple<double, double, double> containing scale, normalized
    * real shift, and normalized real-axis anisotropy.
    */
  constexpr std::tuple<double, double, double> getScaleFromBase();

   /** \brief Shifts and anisotropically transforms a base point without scale.
    * @param lp_base LejaPoint in the base-point set.
    * @param scale_params Transform tuple returned by getScaleFromBase().
    * @return Transformed LejaPoint with its LpType preserved.
    */
  constexpr LejaPoint transformLejaPoint(const LejaPoint& lp_base, const std::tuple<const double, const double, const double>& scale_params);

   /** \brief Sets a positive int polynomial degree and rebuilds base points.
    * @param expansionOrder Positive number of interpolation updates.
    */
  void setExpansionOrder(const int expansionOrder);

   /** \brief Sets ellipse parameters used to transform base Leja points.
    * @pre a <= b and c >= 0.
    * @param a double real-axis minimum parameter.
    * @param b double real-axis maximum parameter.
    * @param c double imaginary semi-axis parameter.
    */
  void setLejaEllipse(const double a, const double b, const double c);

   /** \brief Returns {a, b, c} as Teuchos::Tuple<double, 3>. */
  Teuchos::Tuple<double, 3> getLejaEllipse();

   /** \brief Estimates spectral bounds and updates Leja ellipse parameters.
    *
    * The PhiLinearSolver uses Krylov-Schur or its dense fallback to estimate
    * the real extrema and maximum imaginary magnitude.  The configured safety
    * factor is then applied to a and c only.
    */
  void adaptEvaluator() override;

   /** \brief Selects the int divided-difference algorithm identifier.
    * @param ddMethod 0 recurrence, 1 complex Taylor, 2 dd_phi, or 3 real
    * Taylor; values other than 0--2 select the real Taylor implementation.
    */
  void setDividedDifferenceMethod(const int ddMethod);

   /** \brief Returns one shifted and scaled base Leja point.
    * @param lp_idx Nonnegative index into the stored base-point array.
    * @return LejaPoint transformed to the current ellipse and multiplied by
    * its scalar scale.
    */
  LejaPoint getLpSc(const int lp_idx);

   /** \brief Computes Newton divided differences for the selected method.
    * @param phi_order int phi-function index; the direct evaluation path uses
    * 0.
    * @param cdt Scalar fractional time step used to scale the ellipse.
    * @param lejaOrder int number of coefficients to return.
    * @return Teuchos::ArrayRCP<Scalar> of length @p lejaOrder.
    */
  Teuchos::ArrayRCP<Scalar> getDividedDiffs(const int phi_order, const Scalar cdt, const int lejaOrder);

   /** \brief Returns the configured int polynomial degree. */
  int getExpansionOrder() const { return expansionOrder_; }

  protected:
   /** \brief Applies the Leja polynomial for the matrix exponential in place.
    *
    * @pre phi_order is zero.
    * @param phi_order int phi-function index; only 0 is accepted.
    * @param L Const Teuchos::RCP to the time-scaled Thyra::LinearOpBase<Scalar>.
    * @param v Nonnull Teuchos::Ptr to the input vector, overwritten by the
    * interpolated exponential action.
    * @param cdt Scalar fractional time step used when computing divided
    * differences.
    * @return Thyra::SolveStatus<Scalar> with the final update norm.
    */
   Thyra::SolveStatus<Scalar> computeLinOpPhi(
       const int phi_order,
       const Teuchos::RCP<const Thyra::LinearOpBase<Scalar>> L,
       const Teuchos::Ptr<Thyra::VectorBase<Scalar>> v,
       const Scalar cdt=1.0
     ) override;

private:
  /** \brief int polynomial degree; base storage is rebuilt to support it. */
  int expansionOrder_;
  /** \brief int divided-difference method selected by "Leja DD Method". */
  int ddMethod_;

  /** \brief double safety factor applied to adapted a and c bounds. */
  double leja_sf_;
  /** \brief double absolute convergence threshold for the polynomial update. */
  double leja_tol_;
  /** \brief double minimum real ellipse bound. */
  double leja_a_;
  /** \brief double maximum real ellipse bound. */
  double leja_b_;
  /** \brief double maximum imaginary-magnitude ellipse bound. */
  double leja_c_;

#ifdef TEMPUS_TEUCHOS_TIME_MONITOR
  Teuchos::RCP<Teuchos::Time> timerDD_, timerLinOp_, timerPhi_;
#endif

  // Compile-time indicator that Scalar is a complex type.
  //static const bool isComplex = Teuchos::ScalarTraits<Scalar>::isComplex;

   /** \brief Initializes stored base points.
    * @param lejaOrder Requested storage length, or a negative value to derive
    * it from expansionOrder_.
    */
  void initLejaPointsBase(const int lejaOrder = -1);

   /** \brief Computes coefficients by complex Taylor scaling-and-squaring. */
  Teuchos::ArrayRCP<Scalar> getDividedDiffsTS(const int phi_order, const Scalar cdt, const int lejaOrder);

   /** \brief Computes coefficients by real Taylor scaling-and-squaring. */
  Teuchos::ArrayRCP<Scalar> getDividedDiffsTSR(const int phi_order, const Scalar cdt, const int lejaOrder);

   /** \brief Computes coefficients by recurrence. */
  Teuchos::ArrayRCP<Scalar> getDividedDiffsRC(const int phi_order, const Scalar cdt, const int lejaOrder);

   /** \brief Computes coefficients using the Zivcovich (2019) dd_phi method.
    *
    * Ref: F. Zivcovich, "Fast and accurate computation of divided differences
    * for analytic functions, with an application to the exponential function,"
    * Dolomites Research Notes on Approximation, 12, 2019.
    */
  Teuchos::ArrayRCP<Scalar> getDividedDiffsPhi(const int phi_order, const Scalar cdt, const int lejaOrder);

   /** \brief Teuchos::ArrayRCP of stored base points for the polynomial. */
  Teuchos::ArrayRCP<LejaPoint> lejaPointsBase_;
};

/** \brief Creates a Leja evaluator and applies an optional parameter list.
 *
 * @tparam Scalar Scalar type of the evaluator.
 * @param pList Teuchos::RCP to configuration, or Teuchos::null.
 * @return Teuchos::RCP owning the new PhiEvaluatorLeja<Scalar>.
 */
template <class Scalar>
Teuchos::RCP<PhiEvaluatorLeja<Scalar> > createPhiEvaluatorLeja(
    Teuchos::RCP<Teuchos::ParameterList> pList);

}  // namespace Tempus

#endif  // Tempus_PhiEvaluatorLeja_decl_hpp
