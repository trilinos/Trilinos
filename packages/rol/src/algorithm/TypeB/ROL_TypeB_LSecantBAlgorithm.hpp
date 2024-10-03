// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_TYPEB_LSECANTBALGORITHM_HPP
#define ROL_TYPEB_LSECANTBALGORITHM_HPP

#include "ROL_TypeB_Algorithm.hpp"
#include "ROL_SecantFactory.hpp"
#include "ROL_NullSpaceOperator.hpp"
#include "ROL_ReducedLinearConstraint.hpp"

/** \class ROL::TypeB::LSecantBAlgorithm
    \brief Provides an interface to run the line-search algorithm of Byrd, Lu,
           Nocedal and Zhu (similar to L-BFGS-B).
*/

namespace ROL {
namespace TypeB {

template<typename Real>
class LSecantBAlgorithm : public TypeB::Algorithm<Real> {
private:
  // ITERATION FLAGS/INFORMATION
  int SPflag_;              ///< Subproblem solver termination flag
  int SPiter_;              ///< Subproblem solver iteration count

  // SECANT INFORMATION
  ESecant esec_;             ///< Secant type (default: Limited-Memory BFGS)
  bool useSecantPrecond_;    ///< Flag to use secant as a preconditioner (default: false)
  bool useSecantHessVec_;    ///< Flag to use secant as Hessian (default: false)
  Ptr<Secant<Real>> secant_; ///< Secant object

  // TRUNCATED CG INFORMATION
  Real tol1_; ///< Absolute tolerance for truncated CG (default: 1e-4)
  Real tol2_; ///< Relative tolerance for truncated CG (default: 1e-2)
  int maxit_; ///< Maximum number of CG iterations (default: 20)

  // ALGORITHM SPECIFIC PARAMETERS
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

  unsigned verbosity_; ///< Output level (default: 0)
  bool writeHeader_;   ///< Flag to write header at every iteration

  bool hasEcon_;                            ///< Flag signifies if equality constraints exist
  Ptr<ReducedLinearConstraint<Real>> rcon_; ///< Equality constraint restricted to current active variables
  Ptr<NullSpaceOperator<Real>> ns_;         ///< Null space projection onto reduced equality constraint Jacobian

  using TypeB::Algorithm<Real>::state_;
  using TypeB::Algorithm<Real>::status_;
  using TypeB::Algorithm<Real>::proj_;

public:
  LSecantBAlgorithm(ParameterList &list, const Ptr<Secant<Real>> &secant = nullPtr);

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
  Real dgpstep(Vector<Real> &s, const Vector<Real> &w,
         const Vector<Real> &x, const Real alpha,
         std::ostream &outStream = std::cout) const;

  // Compute Cauchy point, i.e., the minimizer of q(P(x - alpha*g)-x)
  // subject to the trust region constraint ||P(x - alpha*g)-x|| <= del
  //   s      -- The Cauchy step upon return: Primal optimization space vector
  //   alpha  -- The step length for the Cauchy point upon return
  //   x      -- The anchor vector x (unchanged): Primal optimization space vector
  //   g      -- The (dual) gradient vector g (unchanged): Primal optimization space vector
  //   secant -- Secant Hessian approximation
  //   dwa    -- Dual working array, stores Hessian applied to step
  //   dwa1   -- Dual working array
  Real dcauchy(Vector<Real> &s, Real &alpha, Real &q,
               const Vector<Real> &x, const Vector<Real> &g,
               Secant<Real> &secant,
               Vector<Real> &dwa, Vector<Real> &dwa1,
               std::ostream &outStream = std::cout);

  // Perform line search to determine beta such that
  // f(x + beta*s) <= f(x) _ mu0*beta*g's for mu0 in (0,1)
  //   x      -- The anchor vector x, upon return x = x + beta*s: Primal optimization space vector
  //   s      -- The direction vector s, upon return s = beta*s;  Primal optimization space vector
  //   fnew   -- New objective function value
  //   beta   -- Line search step length
  //   fold   -- Old objective function value
  //   gs     -- Gradient applied to the unscaled step s;
  //   obj    -- Objective function
  //   pwa    -- Primal working array
  Real dsrch(Vector<Real> &x, Vector<Real> &s, Real &fnew, Real &beta,
             Real fold, Real gs, Objective<Real> &obj,
             Vector<Real> &pwa, std::ostream &outStream = std::cout);

  // Solve the quadratic subproblem: minimize q(w) subject to the
  // linear constraints using the Conjugate Gradients algorithm
  //   w       -- The step vector to be computed
  //   iflag   -- Termination flag
  //   iter    -- Number of CG iterations
  //   secant  -- Secant Hessian approximation
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
  Real dpcg(Vector<Real> &w, int &iflag, int &iter,
            const Vector<Real> &g, const Vector<Real> &x,
            Secant<Real> &secant, BoundConstraint<Real> &bnd,
            const Real tol, const Real stol, const int itermax,
            Vector<Real> &p, Vector<Real> &q, Vector<Real> &r,
            Vector<Real> &t, Vector<Real> &pwa, Vector<Real> &dwa,
            std::ostream &outStream = std::cout) const;

  void applyFreeHessian(Vector<Real> &hv,
                       const Vector<Real> &v,
                       const Vector<Real> &x,
                       Secant<Real> &secant,
                       BoundConstraint<Real> &bnd,
                       Real &tol,
                       Vector<Real> &pwa) const;

  void applyFreePrecond(Vector<Real> &hv,
                        const Vector<Real> &v,
                        const Vector<Real> &x,
                        Secant<Real> &secant,
                        BoundConstraint<Real> &bnd,
                        Real &tol,
                        Vector<Real> &dwa,
                        Vector<Real> &pwa) const;

}; // class ROL::TypeB::LSecantBAlgorithm

} // namespace TypeB
} // namespace ROL

#include "ROL_TypeB_LSecantBAlgorithm_Def.hpp"

#endif
