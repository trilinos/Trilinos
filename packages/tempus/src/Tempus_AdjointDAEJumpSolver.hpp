//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Tempus_AdjointDAEJumpSolver_hpp
#define Tempus_AdjointDAEJumpSolver_hpp

#include <algorithm>
#include <string>

#include "Teuchos_Array.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Thyra_LinearOpWithSolveBase.hpp"
#include "Thyra_ModelEvaluator.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"

#include "Tempus_AdjointSensitivityUtils.hpp"
#include "Tempus_DAEStructure.hpp"

namespace Tempus {

/** \brief Solver for the adjoint jumps of an index-1 DAE.
 *
 *  For a "Final Time" or "Sum Over Times" response g, the adjoint of an
 *  index-1 DAE jumps at each response time t_k by the differential
 *  components of the solution z of
 *      C^T z = dg/dx^T,   C = P_a A + P_d E,
 *  with E = df/dx_dot, A = df/dx at the forward state at t_k, and P_a, P_d
 *  the masks of the algebraic and differential equations (see DAEStructure),
 *  while its algebraic components give the term -df/dp^T (P_a z) of the
 *  sensitivity (see IntegratorAdjointSensitivity).
 *
 *  C is not assembled.  The solver uses only evaluations of the adjoint model,
 *  which returns the transposed operators: E^T = W(1,0), A^T = W(0,1), and
 *  W(alpha,1) = (alpha E + A)^T.  With D_alpha = P_a + P_d/alpha,
 *      W(alpha,1) D_alpha = C^T + A^T P_d / alpha,
 *  since the algebraic rows of E are zero.  Hence P = W(alpha,1) D_alpha is a
 *  preconditioner of C^T, exact as alpha -> infinity (C is the limit of the
 *  row-scaled step matrix of backward Euler as the time step goes to zero),
 *  and the solver performs the preconditioned Richardson iterations
 *      z <- z + D_alpha^{-1} W(alpha,1)^{-1} (b - C^T z),
 *  with the linear solver of the adjoint model, starting from
 *  alpha = "Initial Alpha Factor" / dt, where dt is the time step size at
 *  t_k.  The iteration matrix is P^{-1} A^T P_d / alpha: whenever an
 *  iteration reduces the residual by a factor larger than "Minimum Residual
 *  Reduction", alpha is multiplied by "Alpha Growth Factor" (and the solver
 *  of W(alpha,1) is rebuilt).
 *
 *  Options (sublist "Algebraic Jump Solver" of "Sensitivities"):
 *  <ul>
 *    <li> "Relative Tolerance" (default: 1e-10) on max_j |r_j|/|b_j| over
 *         the columns j.
 *    <li> "Maximum Iterations Per Alpha" (default: 20)
 *    <li> "Initial Alpha Factor" (default: 1.0)
 *    <li> "Alpha Growth Factor" (default: 10.0)
 *    <li> "Maximum Alpha Increases" (default: 8)
 *    <li> "Minimum Residual Reduction" (default: 0.5)
 *    <li> "Verbose" (default: false) Print the residual at each iteration.
 *  </ul>
 */
template <class Scalar>
class AdjointDAEJumpSolver {
 public:
  typedef Thyra::VectorBase<Scalar> Vector;
  typedef Thyra::MultiVectorBase<Scalar> MultiVector;
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Mag;

  /// Name of the "Sensitivities" sublist with the options of this solver
  static const std::string& sublistName()
  {
    static const std::string name("Algebraic Jump Solver");
    return name;
  }

  static Teuchos::RCP<const Teuchos::ParameterList> getValidParameters()
  {
    Teuchos::RCP<Teuchos::ParameterList> pl =
        Teuchos::parameterList(sublistName());
    pl->set<double>("Relative Tolerance", 1.0e-10, "Tolerance on max_j |b_j - C^T z_j| / |b_j|");
    pl->set<int>("Maximum Iterations Per Alpha", 20, "Maximum number of Richardson iterations for each alpha");
    pl->set<double>("Initial Alpha Factor", 1.0, "Initial alpha, times the time step size at the response time");
    pl->set<double>("Alpha Growth Factor", 10.0, "Factor by which alpha is increased when the iterations stagnate");
    pl->set<int>("Maximum Alpha Increases", 8, "Maximum number of increases of alpha");
    pl->set<double>("Minimum Residual Reduction", 0.5,
                    "alpha is increased when an iteration reduces the residual by a factor larger than this");
    pl->set<bool>("Verbose", false, "Print the residual at each iteration");
    return pl;
  }

  /// Constructor.  pl (optional) is the "Algebraic Jump Solver" sublist.
  AdjointDAEJumpSolver(
      const Teuchos::RCP<const DAEStructure<Scalar> >& dae,
      const Teuchos::RCP<const Teuchos::ParameterList>& pl = Teuchos::null)
    : dae_(dae)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(
        dae_ == Teuchos::null, std::invalid_argument,
        "Error - Tempus::AdjointDAEJumpSolver: null DAE structure.\n");
    Teuchos::ParameterList params(sublistName());
    if (pl != Teuchos::null) params = *pl;
    params.validateParametersAndSetDefaults(*getValidParameters());
    rel_tol_             = params.get<double>("Relative Tolerance");
    max_iters_           = params.get<int>("Maximum Iterations Per Alpha");
    alpha_factor_        = params.get<double>("Initial Alpha Factor");
    alpha_growth_        = params.get<double>("Alpha Growth Factor");
    max_alpha_increases_ = params.get<int>("Maximum Alpha Increases");
    min_reduction_       = params.get<double>("Minimum Residual Reduction");
    verbose_             = params.get<bool>("Verbose");
    TEUCHOS_TEST_FOR_EXCEPTION(
        alpha_factor_ <= 0.0 || alpha_growth_ <= 1.0 || max_iters_ < 1 ||
            max_alpha_increases_ < 0 || min_reduction_ <= 0.0 ||
            min_reduction_ >= 1.0,
        std::invalid_argument,
        "Error - Tempus::AdjointDAEJumpSolver: invalid options (need "
        "\"Initial Alpha Factor\" > 0, \"Alpha Growth Factor\" > 1, "
        "\"Maximum Iterations Per Alpha\" >= 1, \"Maximum Alpha Increases\" "
        ">= 0, 0 < \"Minimum Residual Reduction\" < 1).\n");
  }

  /** \brief Solve C^T z = b for each column of b.
   *
   *  \param adjoint_model  The adjoint model (W_op returns the transpose)
   *  \param inArgs         Forward state at the response time (x, x_dot, t)
   *  \param dt             Time step size at the response time
   *  \param b              Right-hand side, dg/dx^T (x space)
   *  \param z              Solution (residual space), overwritten
   *
   *  Returns true if the iterations converged.
   */
  bool solve(const Thyra::ModelEvaluator<Scalar>& adjoint_model,
             const Thyra::ModelEvaluatorBase::InArgs<Scalar>& inArgs,
             const Scalar dt, const MultiVector& b, MultiVector& z) const
  {
    using Teuchos::RCP;
    typedef Teuchos::ScalarTraits<Scalar> ST;
    typedef Thyra::ModelEvaluatorBase MEB;

    TEUCHOS_TEST_FOR_EXCEPTION(
        !(dt > Scalar(0.0)), std::logic_error,
        "Error - Tempus::AdjointDAEJumpSolver: the time step size must be "
        "positive (dt = "
            << dt << ").\n");
    RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

    // E^T and A^T at the forward state
    RCP<Thyra::LinearOpBase<Scalar> > ET = adjoint_model.create_W_op();
    RCP<Thyra::LinearOpBase<Scalar> > AT = adjoint_model.create_W_op();
    {
      MEB::InArgs<Scalar> in    = inArgs;
      MEB::OutArgs<Scalar> outs = adjoint_model.createOutArgs();
      in.set_alpha(ST::one());
      in.set_beta(ST::zero());
      outs.set_W_op(ET);
      adjoint_model.evalModel(in, outs);
      in.set_alpha(ST::zero());
      in.set_beta(ST::one());
      outs.set_W_op(AT);
      adjoint_model.evalModel(in, outs);
    }

    const int num_cols = b.domain()->dim();
    Teuchos::Array<Mag> b_norms(num_cols);
    Thyra::norms_2(b, b_norms());
    Thyra::assign(Teuchos::ptrFromRef(z), ST::zero());
    if (*std::max_element(b_norms.begin(), b_norms.end()) == Teuchos::ScalarTraits<Mag>::zero())
      return true;

    RCP<MultiVector> r = b.clone_mv();  // residual, x space
    RCP<MultiVector> w = z.clone_mv();  // correction, residual space
    RCP<Vector> d      = Thyra::createMember(dae_->getAlgebraicMask()->space());

    Scalar alpha   = Scalar(alpha_factor_) / dt;
    Mag res        = residual(*ET, *AT, b, z, b_norms, *r);
    int num_solves = 0;
    for (int level = 0; level <= max_alpha_increases_; ++level) {
      // Solver of W(alpha,1) = (alpha E + A)^T, and D_alpha^{-1} = P_a + alpha P_d
      RCP<Thyra::LinearOpWithSolveBase<Scalar> > W =
          AdjointSensitivityUtils::createAdjointLinearOpWithSolve(adjoint_model, inArgs, alpha, ST::one());
      Thyra::assign(d.ptr(), *dae_->getAlgebraicMask());
      Thyra::Vp_StV(d.ptr(), alpha, *dae_->getDifferentialMask());

      for (int it = 0; it < max_iters_; ++it) {
        if (res <= Mag(rel_tol_)) {
          if (verbose_)
            *out << "Tempus::AdjointDAEJumpSolver: converged, residual " << res
                 << ", alpha*dt = " << alpha * dt << ", " << num_solves
                 << " solves" << std::endl;
          return true;
        }
        // z <- z + D_alpha^{-1} W(alpha,1)^{-1} r
        Thyra::assign(w.ptr(), ST::zero());
        W->solve(Thyra::NOTRANS, *r, w.ptr());
        ++num_solves;
        for (int j = 0; j < num_cols; ++j)
          Thyra::ele_wise_scale(*d, w->col(j).ptr());
        Thyra::update(ST::one(), *w, Teuchos::ptrFromRef(z));

        const Mag res_old = res;
        res               = residual(*ET, *AT, b, z, b_norms, *r);
        if (verbose_)
          *out << "Tempus::AdjointDAEJumpSolver: alpha*dt = " << alpha * dt
               << ", iteration " << it + 1 << ", residual " << res
               << std::endl;
        if (res > Mag(min_reduction_) * res_old) break;  // increase alpha
      }
      if (res <= Mag(rel_tol_)) return true;
      alpha *= Scalar(alpha_growth_);
    }
    *out << "WARNING - Tempus::AdjointDAEJumpSolver: no convergence, relative "
         << "residual " << res << " > " << rel_tol_ << " after " << num_solves
         << " solves (final alpha*dt = " << alpha * dt / Scalar(alpha_growth_)
         << ")." << std::endl;
    return false;
  }

 private:
  /// r = b - C^T z, C^T z = A^T (P_a z) + E^T (P_d z); returns
  /// max_j |r_j|/|b_j| over the columns with b_j != 0
  Mag residual(const Thyra::LinearOpBase<Scalar>& ET,
               const Thyra::LinearOpBase<Scalar>& AT, const MultiVector& b,
               const MultiVector& z, const Teuchos::Array<Mag>& b_norms,
               MultiVector& r) const
  {
    typedef Teuchos::ScalarTraits<Scalar> ST;
    Teuchos::RCP<MultiVector> za = z.clone_mv();
    dae_->zeroDifferentialComponents(*za);
    Teuchos::RCP<MultiVector> zd = z.clone_mv();
    dae_->zeroAlgebraicComponents(*zd);
    Thyra::assign(Teuchos::ptrFromRef(r), b);
    AT.apply(Thyra::NOTRANS, *za, Teuchos::ptrFromRef(r), -ST::one(), ST::one());
    ET.apply(Thyra::NOTRANS, *zd, Teuchos::ptrFromRef(r), -ST::one(), ST::one());

    Teuchos::Array<Mag> r_norms(b_norms.size());
    Thyra::norms_2(r, r_norms());
    Mag res = Teuchos::ScalarTraits<Mag>::zero();
    for (int j = 0; j < static_cast<int>(b_norms.size()); ++j)
      if (b_norms[j] > Teuchos::ScalarTraits<Mag>::zero())
        res = std::max(res, r_norms[j] / b_norms[j]);
    return res;
  }

  Teuchos::RCP<const DAEStructure<Scalar> > dae_;
  double rel_tol_;
  int max_iters_;
  double alpha_factor_;
  double alpha_growth_;
  int max_alpha_increases_;
  double min_reduction_;
  bool verbose_;
};

}  // namespace Tempus

#endif  // Tempus_AdjointDAEJumpSolver_hpp
