// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEB_KELLEYSACHSALGORITHM_HPP
#define ROL_TYPEB_KELLEYSACHSALGORITHM_HPP

#include "ROL_TypeB_Algorithm.hpp"
#include "ROL_TrustRegionModel_U.hpp"
#include "ROL_TrustRegionUtilities.hpp"
#include "ROL_NullSpaceOperator.hpp"
#include "ROL_ReducedLinearConstraint.hpp"

/** \class ROL::TypeB::KelleySachsAlgorithm
    \brief Provides an interface to run the trust-region algorithm of Kelley and Sachs.
*/

namespace ROL {
namespace TypeB {

template<typename Real>
class KelleySachsAlgorithm : public TypeB::Algorithm<Real> {
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

  // ITERATION FLAGS/INFORMATION
  TRUtils::ETRFlag TRflag_; ///< Trust-region exit flag
  int SPflag_;              ///< Subproblem solver termination flag
  int SPiter_;              ///< Subproblem solver iteration count

  // SECANT INFORMATION
  ESecant esec_;          ///< Secant type (default: Limited-Memory BFGS)
  bool useSecantPrecond_; ///< Flag to use secant as a preconditioner (default: false)
  bool useSecantHessVec_; ///< Flag to use secant as Hessian (default: false)

  // TRUNCATED CG INFORMATION
  int maxit_; ///< Maximum number of CG iterations (default: 20)
  Real tol1_; ///< Absolute tolerance for truncated CG (default: 1e-4)
  Real tol2_; ///< Relative tolerance for truncated CG (default: 1e-2)

  // ALGORITHM SPECIFIC PARAMETERS
  int minit_;   ///< Maximum number of minor (subproblem solve) iterations (default: 10)
  Real mu0_;    ///< Sufficient decrease parameter (default: 1e-2)
  Real mu1_;    ///< Sufficient decrease parameter postsmoothing (default: 0.9999)
  Real eps0_;   ///< Epsilon binding set tolerance (default: 1e-3)
  Real beta_;   ///< Projected search backtracking rate (default: 1e-2)
  Real alpha0_; ///< Initial step size for projected search (default: 1)

  mutable int nhess_;  ///< Number of Hessian applications
  unsigned verbosity_; ///< Output level (default: 0)
  bool writeHeader_;   ///< Flag to write header at every iteration

  using TypeB::Algorithm<Real>::state_;
  using TypeB::Algorithm<Real>::status_;
  using TypeB::Algorithm<Real>::proj_;

public:
  KelleySachsAlgorithm(ParameterList &list, const Ptr<Secant<Real>> &secant = nullPtr);

  using TypeB::Algorithm<Real>::run;
  void run( Vector<Real>          &x,
            const Vector<Real>    &g, 
            Objective<Real>       &obj,
            BoundConstraint<Real> &bnd,
            std::ostream          &outStream = std::cout) override;

  void run( Problem<Real> &problem,
            std::ostream  &outStream = std::cout ) override;

  void run( Vector<Real>          &x,
            const Vector<Real>    &g,
            Objective<Real>       &obj,
            BoundConstraint<Real> &bnd,
            Constraint<Real>      &linear_econ,
            Vector<Real>          &linear_emul,
            const Vector<Real>    &linear_eres,
            std::ostream          &outStream = std::cout ) override;

  void writeHeader( std::ostream& os ) const override;

  void writeName( std::ostream& os ) const override;

  void writeOutput( std::ostream& os, const bool write_header = false ) const override;

private:
  void initialize(Vector<Real>          &x,
                  const Vector<Real>    &g,
                  Objective<Real>       &obj,
                  BoundConstraint<Real> &bnd,
                  std::ostream &outStream = std::cout);

  // Compute sigma such that ||x+sigma*p||_inv(M) = del.  This is called
  // if trpcg detects negative curvature or if the step violates
  // the trust region bound
  //   xtx -- The dot product <x, inv(M)x> (unchanged)
  //   ptp -- The dot product <p, inv(M)p> (unchanged)
  //   ptx -- The dot product <p, inv(M)x> (unchanged)
  //   del -- Trust region radius (unchanged)
  Real trqsol(const Real xtx, const Real ptp, const Real ptx, const Real del) const;

  // Solve the trust region subproblem: minimize q(w) subject to the
  // trust region constraint ||w||_inv(M) <= del using the Steihaug-Toint
  // Conjugate Gradients algorithm
  //   w       -- The step vector to be computed
  //   iflag   -- Termination flag
  //   iter    -- Number of CG iterations
  //   g       -- Free components of gradient (unchanged)
  //   x       -- Current iterate (unchanged)
  //   g0      -- Gradient (unchanged)
  //   del     -- Trust region radius (unchanged)
  //   model   -- Trust region model
  //   bnd     -- Bound constraint used to remove active variables
  //   eps     -- Epsilon binding set tolerance
  //   tol     -- Residual stopping tolerance (unchanged)
  //   stol    -- Preconditioned residual stopping tolerance (unchanged)
  //   itermax -- Maximum number of iterations
  //   p       -- Primal working array that stores the CG step
  //   q       -- Dual working array that stores the Hessian applied to p
  //   r       -- Primal working array that stores the preconditioned residual
  //   t       -- Dual working array that stores the residual
  //   pwa     -- Primal working array that stores the pruned vector in hessVec
  //   dwa     -- Dual working array that stores the pruned vector in precond
  Real trpcg(Vector<Real> &w, int &iflag, int &iter,
             const Vector<Real> &g, const Vector<Real> &x,
             const Vector<Real> &g0,
             const Real del, TrustRegionModel_U<Real> &model,
             BoundConstraint<Real> &bnd, Real eps,
             const Real tol, const int itermax,
             Vector<Real> &p, Vector<Real> &q, Vector<Real> &r,
             Vector<Real> &t, Vector<Real> &pwa, Vector<Real> &dwa,
             std::ostream &outStream = std::cout) const;

  void applyFreeHessian(Vector<Real> &hv,
                       const Vector<Real> &v,
                       const Vector<Real> &x,
                       const Vector<Real> &g,
                       Real eps,
                       TrustRegionModel_U<Real> &model,
                       BoundConstraint<Real> &bnd,
                       Real &tol,
                       Vector<Real> &pwa,
                       Vector<Real> &dwa) const;

  void applyFreePrecond(Vector<Real> &hv,
                        const Vector<Real> &v,
                        const Vector<Real> &x,
                         const Vector<Real> &g,
                         Real eps,
                        TrustRegionModel_U<Real> &model,
                        BoundConstraint<Real> &bnd,
                        Real &tol,
                        Vector<Real> &pwa,
                        Vector<Real> &dwa) const;

}; // class ROL::TypeB::KelleySachsAlgorithm

} // namespace TypeB
} // namespace ROL

#include "ROL_TypeB_KelleySachsAlgorithm_Def.hpp"

#endif
