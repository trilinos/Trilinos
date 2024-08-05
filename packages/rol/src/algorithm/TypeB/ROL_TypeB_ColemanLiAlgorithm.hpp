// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEB_COLEMANLIALGORITHM_HPP
#define ROL_TYPEB_COLEMANLIALGORITHM_HPP

#include "ROL_TypeB_Algorithm.hpp"
#include "ROL_TrustRegionModel_U.hpp"
#include "ROL_TrustRegionUtilities.hpp"
#include "ROL_NullSpaceOperator.hpp"
#include "ROL_ReducedLinearConstraint.hpp"

/** \class ROL::TypeB::ColemanLiAlgorithm
    \brief Provides an interface to run the affine-scaling trust-region algorithm of Coleman and Li.
*/

namespace ROL {
namespace TypeB {

template<typename Real>
class ColemanLiAlgorithm : public TypeB::Algorithm<Real> {
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
  Real mu0_;       ///< Sufficient decrease parameter (default: 1e-2)
  Real spexp_;     ///< Relative tolerance exponent for subproblem solve (default: 1, range: [1,2])
  Real alphaMax_;  ///< Maximum value of relaxation parameter (default: 0.999)

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
  ColemanLiAlgorithm(ParameterList &list, const Ptr<Secant<Real>> &secant = nullPtr);

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

  // Compute the projected step s = P(x + alpha*w) - x
  // Returns the norm of the projected step s
  //    s     -- The projected step upon return
  //    w     -- The direction vector w (unchanged)
  //    x     -- The anchor vector x (unchanged)
  //    alpha -- The step size (unchanged)
  //    bnd   -- The bound constraint
  Real dgpstep(Vector<Real> &s, const Vector<Real> &w,
         const Vector<Real> &x, const Real alpha,
         std::ostream &outStream = std::cout) const;

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
  //   g       -- Vector containing gradient (dual space)
  //   x       -- Vector containing iterate (primal space)
  //   gdual   -- Vector containing primal gradient (primal space)
  //   del     -- Trust region radius (unchanged)
  //   model   -- Trust region model
  //   bnd     -- Bound constraint used to remove active variables
  //   tol     -- Residual stopping tolerance (unchanged)
  //   stol    -- Preconditioned residual stopping tolerance (unchanged)
  //   p       -- Primal working array that stores the CG step
  //   q       -- Dual working array that stores the Hessian applied to p
  //   r       -- Primal working array that stores the preconditioned residual
  //   t       -- Dual working array that stores the residual
  //   pwa1    -- Primal working array that stores the pruned vector in hessVec
  //   pwa2    -- Primal working array that stores the pruned vector in hessVec
  //   dwa     -- Dual working array that stores the pruned vector in precond
  Real dtrpcg(Vector<Real> &w, int &iflag, int &iter,
              const Vector<Real> &g, const Vector<Real> &x, const Vector<Real> &gdual,
              const Real del, TrustRegionModel_U<Real> &model, BoundConstraint<Real> &bnd,
              const Real tol, const Real stol,
              Vector<Real> &p, Vector<Real> &q, Vector<Real> &r,
              Vector<Real> &t, Vector<Real> &pwa1, Vector<Real> &pwa2,
              Vector<Real> &dwa,
              std::ostream &outStream = std::cout) const;

  void applyC(Vector<Real> &Cv,
             const Vector<Real> &v,
             const Vector<Real> &x,
             const Vector<Real> &g,
             BoundConstraint<Real> &bnd,
             Vector<Real> &pwa) const;

  void applyHessian(Vector<Real> &hv,
                   const Vector<Real> &v,
                   const Vector<Real> &x,
                   const Vector<Real> &g,
                   TrustRegionModel_U<Real> &model,
                   BoundConstraint<Real> &bnd,
                   Real &tol,
                   Vector<Real> &pwa1,
                   Vector<Real> &pwa2) const;

  void applyPrecond(Vector<Real> &hv,
                    const Vector<Real> &v,
                    const Vector<Real> &x,
                    const Vector<Real> &g,
                    TrustRegionModel_U<Real> &model,
                    BoundConstraint<Real> &bnd,
                    Real &tol,
                    Vector<Real> &dwa,
                    Vector<Real> &pwa) const;

}; // class ROL::TypeB::ColemanLiAlgorithm

} // namespace TypeB
} // namespace ROL

#include "ROL_TypeB_ColemanLiAlgorithm_Def.hpp"

#endif
