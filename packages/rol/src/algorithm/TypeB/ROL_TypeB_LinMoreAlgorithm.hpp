// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEB_LINMOREALGORITHM_HPP
#define ROL_TYPEB_LINMOREALGORITHM_HPP

#include "ROL_TypeB_Algorithm.hpp"
#include "ROL_TrustRegionModel_U.hpp"
#include "ROL_TrustRegionUtilities.hpp"
#include "ROL_NullSpaceOperator.hpp"
#include "ROL_ReducedLinearConstraint.hpp"

/** \class ROL::TypeB::LinMoreAlgorithm
    \brief Provides an interface to run the trust-region algorithm of Lin and More.
*/

namespace ROL {
namespace TypeB {

template<typename Real>
class LinMoreAlgorithm : public TypeB::Algorithm<Real> {
private:
  Ptr<TrustRegionModel_U<Real>> model_;  ///< Container for trust-region model

  // TRUST REGION PARAMETERS
  Real delMax_; ///< Maximum trust-region radius (default: ROL_INF)
  Real eta0_;   ///< Step acceptance threshold (default: 0.05)
  Real eta1_;   ///< Radius decrease threshold (default: 0.05)
  Real eta2_;   ///< Radius increase threshold (default: 0.9)
  Real gamma0_; ///< Radius decrease rate (negative rho) (default: 0.0625)
  Real gamma1_; ///< Radius decrease rate (positive rho) (default: 0.25)
  Real gamma2_; ///< Radius increase rate (default: 2.5)
  Real TRsafe_; ///< Safeguard size for numerically evaluating ratio (default: 1e2)
  Real eps_;    ///< Safeguard for numerically evaluating ratio
  bool interpRad_; ///< Interpolate the trust-region radius if ratio is negative (default: false)

  // NONMONOTONE PARAMETER
  bool useNM_;
  int storageNM_;

  // ITERATION FLAGS/INFORMATION
  TRUtils::ETRFlag TRflag_; ///< Trust-region exit flag
  int SPflag_;              ///< Subproblem solver termination flag
  int SPiter_;              ///< Subproblem solver iteration count

  // SECANT INFORMATION
  ESecant esec_;          ///< Secant type (default: Limited-Memory BFGS)
  bool useSecantPrecond_; ///< Flag to use secant as a preconditioner (default: false)
  bool useSecantHessVec_; ///< Flag to use secant as Hessian (default: false)

  // TRUNCATED CG INFORMATION
  Real tol1_; ///< Absolute tolerance for truncated CG (default: 1e-4)
  Real tol2_; ///< Relative tolerance for truncated CG (default: 1e-2)
  int maxit_; ///< Maximum number of CG iterations (default: 20)

  // ALGORITHM SPECIFIC PARAMETERS
  int minit_;      ///< Maximum number of minor (subproblem solve) iterations (default: 10)
  Real mu0_;       ///< Sufficient decrease parameter (default: 1e-2)
  Real spexp_;     ///< Relative tolerance exponent for subproblem solve (default: 1, range: [1,2])
  int  redlim_;    ///< Maximum number of Cauchy point reduction steps (default: 10)
  int  explim_;    ///< Maximum number of Cauchy point expansion steps (default: 10)
  Real alpha_;     ///< Initial Cauchy point step length (default: 1.0)
  bool normAlpha_; ///< Normalize initial Cauchy point step length (default: false)
  Real interpf_;   ///< Backtracking rate for Cauchy point computation (default: 1e-1)
  Real extrapf_;   ///< Extrapolation rate for Cauchy point computation (default: 1e1)
  Real qtol_;      ///< Relative tolerance for computed decrease in Cauchy point computation (default: 1-8)
  Real interpfPS_; ///< Backtracking rate for projected search (default: 0.5)
  int pslim_;      ///< Maximum number of projected search steps (default: 20)

  // Inexactness Parameters
  std::vector<bool> useInexact_;
  Real scale0_;
  Real scale1_;
  Real scale_;
  Real omega_;
  Real force_;
  int updateIter_;
  Real forceFactor_;
  Real gtol_;

  mutable int nhess_;  ///< Number of Hessian applications
  unsigned verbosity_; ///< Output level (default: 0)
  bool writeHeader_;   ///< Flag to write header at every iteration

  bool hasEcon_;                            ///< Flag signifies if equality constraints exist
  Ptr<ReducedLinearConstraint<Real>> rcon_; ///< Equality constraint restricted to current active variables
  Ptr<NullSpaceOperator<Real>> ns_;         ///< Null space projection onto reduced equality constraint Jacobian

  using TypeB::Algorithm<Real>::state_;
  using TypeB::Algorithm<Real>::status_;
  using TypeB::Algorithm<Real>::proj_;

public:
  LinMoreAlgorithm(ParameterList &list, const Ptr<Secant<Real>> &secant = nullPtr);

  using TypeB::Algorithm<Real>::run;
  void run( Vector<Real>          &x,
            const Vector<Real>    &g, 
            Objective<Real>       &obj,
            BoundConstraint<Real> &bnd,
            std::ostream          &outStream = std::cout) override;

  void writeHeader( std::ostream& os ) const override;

  void writeName( std::ostream& os ) const override;

  void writeOutput( std::ostream& os, const bool write_header = false ) const override;

private:
  void initialize(Vector<Real>          &x,
                  const Vector<Real>    &g,
                  Objective<Real>       &obj,
                  BoundConstraint<Real> &bnd,
                  std::ostream &outStream = std::cout);

  Real computeValue(Real inTol,
                    Real &outTol,
                    Real pRed,
                    Real &fold,
                    int iter,
                    const Vector<Real> &x,
                    const Vector<Real> &xold,
                    Objective<Real> &obj);

  void computeGradient(const Vector<Real> &x,
                       Vector<Real> &g,
                       Vector<Real> &pwa,
                       Real del,
                       Objective<Real> &obj,
                       bool accept,
                       Real &gtol,
                       Real &gnorm,
                       std::ostream &outStream = std::cout) const;

  // Compute the projected step s = P(x + alpha*w) - x
  // Returns the norm of the projected step s
  //    s     -- The projected step upon return
  //    w     -- The direction vector w (unchanged)
  //    x     -- The anchor vector x (unchanged)
  //    alpha -- The step size (unchanged)
  Real dgpstep(Vector<Real> &s, const Vector<Real> &w,
         const Vector<Real> &x, const Real alpha,
         std::ostream &outStream = std::cout) const;

  // Compute Cauchy point, i.e., the minimizer of q(P(x - alpha*g)-x)
  // subject to the trust region constraint ||P(x - alpha*g)-x|| <= del
  //   s     -- The Cauchy step upon return: Primal optimization space vector
  //   alpha -- The step length for the Cauchy point upon return
  //   x     -- The anchor vector x (unchanged): Primal optimization space vector
  //   g     -- The (dual) gradient vector g (unchanged): Primal optimization space vector
  //   del   -- The trust region radius (unchanged)
  //   model -- Trust region model
  //   dwa   -- Dual working array, stores Hessian applied to step
  //   dwa1  -- Dual working array
  Real dcauchy(Vector<Real> &s, Real &alpha, Real &q,
               const Vector<Real> &x, const Vector<Real> &g,
               const Real del, TrustRegionModel_U<Real> &model,
               Vector<Real> &dwa, Vector<Real> &dwa1,
               std::ostream &outStream = std::cout);

  // Perform projected search to determine beta such that
  // q(P(x + beta*s)-x) <= mu0*g'(P(x + beta*s)-x) for mu0 in (0,1)
  //   x     -- The anchor vector x, upon return x = P(x + beta*s): Primal optimization space vector
  //   s     -- The direction vector s, upon return s = P(x + beta*s) - x: Primal optimization space vector
  //   g     -- The free components of the gradient vector g (unchanged): Primal optimization space vector
  //   model -- Contains objective and bound constraint information
  //   pwa   -- Primal working array
  //   dwa   -- Dual working array
  Real dprsrch(Vector<Real> &x, Vector<Real> &s, Real &q,
               const Vector<Real> &g, TrustRegionModel_U<Real> &model,
               BoundConstraint<Real> &bnd,
               Vector<Real> &pwa, Vector<Real> &dwa,
               std::ostream &outStream = std::cout);

  // Compute sigma such that ||x+sigma*p||_inv(M) = del.  This is called
  // if dtrpcg detects negative curvature or if the step violates
  // the trust region bound
  //   xtx -- The dot product <x, inv(M)x> (unchanged)
  //   ptp -- The dot product <p, inv(M)p> (unchanged)
  //   ptx -- The dot product <p, inv(M)x> (unchanged)
  //   del -- Trust region radius (unchanged)
  Real dtrqsol(const Real xtx, const Real ptp, const Real ptx, const Real del) const;

  // Solve the trust region subproblem: minimize q(w) subject to the
  // trust region constraint ||w||_inv(M) <= del using the Steihaug-Toint
  // Conjugate Gradients algorithm
  //   w       -- The step vector to be computed
  //   iflag   -- Termination flag
  //   iter    -- Number of CG iterations
  //   del     -- Trust region radius (unchanged)
  //   model   -- Trust region model
  //   bnd     -- Bound constraint used to remove active variables
  //   tol     -- Residual stopping tolerance (unchanged)
  //   stol    -- Preconditioned residual stopping tolerance (unchanged)
  //   itermax -- Maximum number of iterations
  //   p       -- Primal working array that stores the CG step
  //   q       -- Dual working array that stores the Hessian applied to p
  //   r       -- Primal working array that stores the preconditioned residual
  //   t       -- Dual working array that stores the residual
  //   pwa     -- Primal working array that stores the pruned vector in hessVec
  //   dwa     -- Dual working array that stores the pruned vector in precond
  Real dtrpcg(Vector<Real> &w, int &iflag, int &iter,
              const Vector<Real> &g, const Vector<Real> &x,
              const Real del, TrustRegionModel_U<Real> &model, BoundConstraint<Real> &bnd,
              const Real tol, const Real stol, const int itermax,
              Vector<Real> &p, Vector<Real> &q, Vector<Real> &r,
              Vector<Real> &t, Vector<Real> &pwa, Vector<Real> &dwa,
              std::ostream &outStream = std::cout) const;

  void applyFreeHessian(Vector<Real> &hv,
                       const Vector<Real> &v,
                       const Vector<Real> &x,
                       TrustRegionModel_U<Real> &model,
                       BoundConstraint<Real> &bnd,
                       Real &tol,
                       Vector<Real> &pwa) const;

  void applyFreePrecond(Vector<Real> &hv,
                        const Vector<Real> &v,
                        const Vector<Real> &x,
                        TrustRegionModel_U<Real> &model,
                        BoundConstraint<Real> &bnd,
                        Real &tol,
                        Vector<Real> &dwa,
                        Vector<Real> &pwa) const;

}; // class ROL::TypeB::LinMoreAlgorithm

} // namespace TypeB
} // namespace ROL

#include "ROL_TypeB_LinMoreAlgorithm_Def.hpp"

#endif
