// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_LINMOREALGORITHM_B_H
#define ROL_LINMOREALGORITHM_B_H

#include "ROL_Algorithm_B.hpp"

/** \class ROL::LinMoreAlgorithm_B
    \brief Provides an interface to run the trust-region algorithm of Lin and More.
*/

namespace ROL {

template<typename Real>
class LinMoreAlgorithm_B : public Algorithm_B<Real> {
private:
  Ptr<Vector<Real>> x_, s_, g_;
  Ptr<Vector<Real>> pwa1_, pwa2_, pwa3_, dwa1_, dwa2_, dwa3_;

  Real tol1_, tol2_, alpha_;
  int maxit_;

  unsigned verbosity_;
  bool printHeader_;

  bool hasEcon_;
  Ptr<ReducedConstraint<Real>> rcon_;
  PtrVector<Real>> ran_;

  class ReducedConstraint : Constraint<Real> {
    private:
      const Ptr<Constraint<Real>> con_;
      const Ptr<BoundConstraint<Real>> bnd_;
      Ptr<const Vector<Real>> x_;
      Ptr<Vector<Real>> prim_;
    public:
      ReducedConstraint(const Ptr<Constraint<Real>> &con,
                        const Ptr<BoundConstraint<Real>> &bnd,
                        const Ptr<const Vector<Real>> &x)
        : con_(con), bnd_(bnd), x_(x), prim_(x->clone()) {}

      void setX(const Ptr<const Vector<Real>> &x) {
        x_ = x;
      }

      void value(Vector<Real> &c, const Vector<Real> &x, Real &tol) {
        const Real zero(0);
        prim_->set(x);
        bnd_->pruneActive(prim_,*x_,zero);
        con_->value(c,*prim_,zero);
      }

      void applyJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
        const Real zero(0);
        prim_->set(v);
        bnd_->pruneActive(prim_,*x_,zero);
        con_->applyJacobian(jv,*prim_,x,zero);
      }

      void applyAdjointJacobian(Vector<Real> &jv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
        const Real zero(0);
        con_->applyAdjointJacobian(jv,v,x,tol);
        bnd_->pruneActive(jv,*x_,zero);
      }

      void applyAdjointHessian(Vector<Real> &ahuv, const Vector<Real> &u, const Vector<Real> &v,
                               const Vector<Real> &x, Real &tol) {
        ahuv.zero();
      }
  };

  using Algorithm_B<Real>::state_;
  using Algorithm_B<Real>::status_;
  using Algorithm_B<Real>::proj_;

public:
  LinMoreAlgorithm_B(ParameterList &list);

  std::vector<std::string> run( Vector<Real>          &x,
                                const Vector<Real>    &g, 
                                Objective<Real>       &obj,
                                BoundConstraint<Real> &bnd,
                                std::ostream          &outStream = std::cout);

  void solve( Vector<Real>           &s,
              Real                   &snorm,
              int                    &iflag,
              int                    &iter,
              const Real              del,
              TrustRegionModel<Real> &model );

  // Compute the projected step s = P(x + alpha*w) - x
  // Returns the norm of the projected step s
  //    s     -- The projected step upon return
  //    w     -- The direction vector w (unchanged)
  //    x     -- The anchor vector x (unchanged)
  //    alpha -- The step size (unchanged)
  Real dgpstep(Vector<Real> &s, const Vector<Real> &w,
         const Vector<Real> &x, const Real alpha) const;

  // Compute Cauchy point, i.e., the minimizer of q(P(x - alpha*g)-x)
  // subject to the trust region constraint ||P(x - alpha*g)-x|| <= del
  //   s     -- The Cauchy point upon return
  //   alpha -- The step length for the Cauchy point upon return
  //   x     -- The anchor vector x (unchanged)
  //   g     -- The (dual) gradient vector g (unchanged)
  //   del   -- The trust region radius (unchanged)
  //   obj   -- Objective function (used for hessVec)
  //   dwa   -- Dual working array
  Real dcauchy(Vector<Real> &s, Real &alpha,
               const Vector<Real> &x, const Vector<Real> &g,
               const Real del, Objective<Real> &obj,
               Vector<Real> &pwa1, Vector<Real> &pwa2, Vector<Real> &dwa,
               std::ostream &outStream = std::cout);

  // Perform projected search to determine beta such that
  // q(P(x + beta*s)-x) <= mu0*g'(P(x + beta*s)-x) for mu0 in (0,1)
  //   x     -- The anchor vector x, upon return x = P(x + beta*s)
  //   s     -- The direction vector s, upon return s = P(x + beta*s) - x
  //   g     -- The free components of the gradient vector g (unchanged)
  //   model -- Contains objective and bound constraint information
  //   pwa   -- Primal working array
  //   dwa   -- Dual working array
  Real dprsrch(Vector<Real> &x, Vector<Real> &s,
               const Vector<Real> &g, Objective<Real> &obj,
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
  //   obj     -- Objective function (used for hessVec)
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
              const Real del, Objective<Real> &obj, BoundConstraint<Real> &bnd,
              const Real tol, const Real stol, const int itermax,
              Vector<Real> &p, Vector<Real> &q, Vector<Real> &r,
              Vector<Real> &t, Vector<Real> &pwa, const Vector<Real> &dwa) const;

  void applyFreeHessian(Vector<Real> &hv,
                       const Vector<Real> &v,
                       const Vector<Real> &x,
                       Objective<Real> &obj,
                       BoundConstraint<Real> &bnd,
                       Real &tol,
                       Vector<Real> &pwa);

  void applyFreePrecond(Vector<Real> &hv,
                        const Vector<Real> &v,
                        const Vector<Real> &x,
                        Objective<Real> &obj,
                        BoundConstraint<Real> &bnd,
                        Real &tol,
                        Vector<Real> &dwa) {

  std::string printHeader( void ) const;

  std::string printName( void ) const;

  std::string print( const bool print_header = false ) const;

}; // class ROL::LinMoreAlgorithm_B

} // namespace ROL

#include "ROL_LinMoreAlgorithm_B_Def.hpp"

#endif
