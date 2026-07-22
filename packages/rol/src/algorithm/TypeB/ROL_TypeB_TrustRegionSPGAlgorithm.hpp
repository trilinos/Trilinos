// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEB_TRUSTREGIONSPGALGORITHM_HPP
#define ROL_TYPEB_TRUSTREGIONSPGALGORITHM_HPP

#include "ROL_TypeB_Algorithm.hpp"
#include "ROL_TrustRegionModel_U.hpp"
#include "ROL_TrustRegionUtilities.hpp"

/** \class ROL::TypeB::TrustRegionSPGAlgorithm
    \brief Provides an interface to run the trust-region algorithm.
*/

namespace ROL {
namespace TypeB {

template<typename Real>
class TrustRegionSPGAlgorithm : public TypeB::Algorithm<Real> {
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
  int maxit_; ///< Maximum number of CG iterations (default: 25)

  // ALGORITHM SPECIFIC PARAMETERS
  bool useNM_;
  int  storageNM_;
  Real mu0_;       ///< Sufficient decrease parameter (default: 1e-2)
  Real spexp_;     ///< Relative tolerance exponent for subproblem solve (default: 1, range: [1,2])
  int  redlim_;    ///< Maximum number of Cauchy point reduction steps (default: 10)
  int  explim_;    ///< Maximum number of Cauchy point expansion steps (default: 10)
  Real alpha_;     ///< Initial Cauchy point step length (default: 1.0)
  bool normAlpha_; ///< Normalize initial Cauchy point step length (default: false)
  Real interpf_;   ///< Backtracking rate for Cauchy point computation (default: 1e-1)
  Real extrapf_;   ///< Extrapolation rate for Cauchy point computation (default: 1e1)
  Real qtol_;      ///< Relative tolerance for computed decrease in Cauchy point computation (default: 1-8)
  Real lambdaMin_;
  Real lambdaMax_;
  Real gamma_;
  int maxSize_;
  bool useMin_;
  bool useNMSP_;
  bool useSimpleSPG_;

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

  using TypeB::Algorithm<Real>::state_;
  using TypeB::Algorithm<Real>::status_;
  using TypeB::Algorithm<Real>::proj_;

public:
  TrustRegionSPGAlgorithm(ParameterList &list, const Ptr<Secant<Real>> &secant = nullPtr);

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
                  Real                   ftol,
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

  void dpsg_simple(Vector<Real> &y, Real &q, Vector<Real> &gmod, const Vector<Real> &x,
                   Real del, TrustRegionModel_U<Real> &model, Vector<Real> &pwa,
                   Vector<Real> &pwa1, Vector<Real> &dwa, std::ostream &outStream = std::cout);

  void dpsg(Vector<Real> &y, Real &q, Vector<Real> &gmod, const Vector<Real> &x,
            Real del, TrustRegionModel_U<Real> &model, Vector<Real> &ymin,
            Vector<Real> &pwa, Vector<Real> &pwa1, Vector<Real> &pwa2,
            Vector<Real> &pwa3, Vector<Real> &pwa4, Vector<Real> &pwa5,
            Vector<Real> &dwa, std::ostream &outStream = std::cout);

  void dproj(Vector<Real> &x, const Vector<Real> &x0, Real del,
            Vector<Real> &y0, Vector<Real> &y1, Vector<Real> &yc,
            Vector<Real> &pwa, std::ostream &outStream = std::cout) const;

}; // class ROL::TypeB::TrustRegionSPGAlgorithm

} // namespace TypeB
} // namespace ROL

#include "ROL_TypeB_TrustRegionSPGAlgorithm_Def.hpp"

#endif
