// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_LINMORE_H
#define ROL_LINMORE_H

/** \class ROL::LinMore
    \brief Provides interface for truncated CG trust-region subproblem solver.
*/

#include "ROL_TrustRegion.hpp"
#include "ROL_LinMoreModel.hpp"
#include "ROL_Elementwise_Function.hpp"
#include "ROL_Elementwise_Reduce.hpp"

namespace ROL { 

template<class Real>
class LinMore : public TrustRegion<Real> {
private:

  Ptr<Vector<Real>> x_, s_, g_;
  Ptr<Vector<Real>> pwa1_, pwa2_, dwa1_, dwa2_;

  Real tol1_, tol2_, alpha_;
  int maxit_;

  unsigned verbosity_;

  class LowerBreakPoint : public Elementwise::BinaryFunction<Real> {
    public:
    Real apply( const Real &x, const Real &y) const {
      const Real zero(0), one(1);
      return (x > zero && y < zero) ? -x/y : -one;
    }
  } lbp_;

  class UpperBreakPoint : public Elementwise::BinaryFunction<Real> {
    public:
    Real apply( const Real &x, const Real &y) const {
      const Real zero(0), one(1);
      return (x > zero && y > zero) ? x/y : -one;
    }
  } ubp_;

  class PositiveMin : public Elementwise::ReductionOp<Real> {
    public:
    void reduce(const Real &input, Real &output) const {
      const Real zero(0);
      output = (input<output && input>zero) ? input : output;
    }
    void reduce( const volatile Real &input, Real volatile &output ) const {
      const Real zero(0);
      output = (input<output && input>zero) ? input : output;
    }
    Real initialValue() const {
      return ROL_INF<Real>();
    }
    Elementwise::EReductionType reductionType() const {
      return Elementwise::REDUCE_MIN;
    }
  } pmin_;

  class PositiveMax : public Elementwise::ReductionOp<Real> {
    public:
    void reduce(const Real &input, Real &output) const {
      const Real zero(0);
      output = (input>output && input>zero) ? input : output;
    }
    void reduce( const volatile Real &input, Real volatile &output ) const {
      const Real zero(0);
      output = (input>output && input>zero) ? input : output;
    }
    Real initialValue() const {
      return static_cast<Real>(0);
    }
    Elementwise::EReductionType reductionType() const {
      return Elementwise::REDUCE_MAX;
    }
  } pmax_;

public:

  // Constructor
  LinMore( ROL::ParameterList &parlist ) : TrustRegion<Real>(parlist), alpha_(1) {
    // Unravel Parameter List
    Real em4(1e-4), em2(1e-2);
    maxit_   = parlist.sublist("General").sublist("Krylov").get("Iteration Limit",20);
    tol1_    = parlist.sublist("General").sublist("Krylov").get("Absolute Tolerance",em4);
    tol2_    = parlist.sublist("General").sublist("Krylov").get("Relative Tolerance",em2);
    // Get verbosity level
    verbosity_ = parlist.sublist("General").get("Print Verbosity", 0);
  }

  void initialize( const Vector<Real> &x, const Vector<Real> &s, const Vector<Real> &g) {
    TrustRegion<Real>::initialize(x,s,g);
    x_ = x.clone(); s_ = x.clone(); g_ = g.clone();
    pwa1_ = x.clone(); pwa2_ = x.clone();
    dwa1_ = g.clone(); dwa2_ = g.clone();
  }

  void run( Vector<Real>           &s,
            Real                   &snorm,
            int                    &iflag,
            int                    &iter,
            const Real              del,
            TrustRegionModel<Real> &model ) {
    const Real zero(0), half(0.5), one(1);
    Real tol0 = std::sqrt(ROL_EPSILON<Real>());
    Real gfnorm(0), gfnormf(0), tol(0), stol(0);
    int dim = s.dimension();
    // Compute Cauchy point (TRON notation: x_ = x[1])
    snorm = dcauchy(*s_,alpha_,*model.getIterate(),model.getGradient()->dual(),
                    del,model,*pwa1_,*pwa2_,*dwa1_); // Solve 1D optimization problem for alpha_
    x_->set(*model.getIterate());                    // TRON notation: model.getIterate() = x[0]
    x_->plus(*s_);                                   // Set x_ = x[0] + alpha_*g
    model.getBoundConstraint()->project(*x_);        // Project x_ onto bounds

    // Model gradient at s = x[1] - x[0]
    s.set(*x_); s.axpy(-one,*model.getIterate()); // s_ = x[i+1]-x[0]
    dynamic_cast<LinMoreModel<Real>&>(model).applyFullHessian(*g_,s,tol0);
    g_->plus(*model.getGradient());
    model.getBoundConstraint()->pruneActive(*g_,*x_,zero);
    gfnorm = g_->norm();
    if (verbosity_ > 0) {
      std::cout << std::endl;
      std::cout << "  Computation of Cauchy point"          << std::endl;
      std::cout << "    Step length (alpha):              " << alpha_     << std::endl;
      std::cout << "    Step length (alpha*g):            " << snorm      << std::endl;
      std::cout << "    Norm of Cauchy point:             " << x_->norm() << std::endl;
      std::cout << "    Norm of free gradient components: " << gfnorm     << std::endl;
    }

    // Main step computation loop
    // There are at most dim iterations since at least one
    // face becomes active at each iteration
    iter = 0;
    for (int i = 0; i < dim; ++i) {
      // Run Truncated CG
      int flagCG = 0, iterCG = 0;
      tol  = tol2_*gfnorm;
      stol = zero;
      snorm = dtrpcg(*s_,flagCG,iterCG,*g_,*x_,del,model,
                     tol,stol,maxit_,
                     *pwa1_,*dwa1_,*pwa2_,*dwa2_);
      iter += iterCG;
      if (verbosity_ > 0) {
        std::cout << std::endl;
        std::cout << "  Computation of CG step"               << std::endl;
        std::cout << "    Number of faces:                  " << dim        << std::endl;
        std::cout << "    Current face (i):                 " << i          << std::endl;
        std::cout << "    CG step length:                   " << snorm      << std::endl;
        std::cout << "    Number of CG iterations:          " << iterCG     << std::endl;
        std::cout << "    CG flag:                          " << flagCG     << std::endl;
        std::cout << "    Total number of iterations:       " << iter       << std::endl;
      }

      // Projected search
      snorm = dprsrch(*x_,*s_,g_->dual(),model,*pwa1_,*dwa1_);
      if (verbosity_ > 0) {
        std::cout << "    Step length (beta*s):             " << snorm      << std::endl;
        std::cout << "    Iterate length:                   " << x_->norm() << std::endl;
      }

      // Model gradient at s = x[i+1] - x[0]
      s.set(*x_); s.axpy(-one,*model.getIterate()); // s_ = x[i+1]-x[0]
      dynamic_cast<LinMoreModel<Real>&>(model).applyFullHessian(*g_,s,tol0);
      g_->plus(*model.getGradient());
      model.getBoundConstraint()->pruneActive(*g_,*x_,zero);
      gfnormf = g_->norm();
      if (verbosity_ > 0) {
        std::cout << std::endl;
        std::cout << "  Update model gradient"                << std::endl;
        std::cout << "    Step length:                      " << s.norm()   << std::endl;
        std::cout << "    Norm of free gradient components: " << gfnormf    << std::endl;
        std::cout << std::endl;
      }

      // Termination check
      if (gfnormf <= tol2_*gfnorm) {
        iflag = 0;
        break;
      }
      else if (iter >= maxit_) {
        iflag = 1;
        break;
      }
      else if (flagCG == 2) {
        iflag = 2;
        break;
      }
      else if (flagCG == 3) {
        iflag = 3;
        break;
      }

      // Update free gradient norm
      gfnorm = gfnormf;
    }
    // Update norm of step and update model predicted reduction
    snorm = s.norm();
    Real gs(0), q(0);
    dynamic_cast<LinMoreModel<Real>&>(model).applyFullHessian(*dwa1_,s,tol0);
    gs = s.dot(model.getGradient()->dual());
    q  = half * s.dot(dwa1_->dual()) + gs;
    TrustRegion<Real>::setPredictedReduction(-q);
  }

private:

  // Compute the projected step s = P(x + alpha*w) - x
  // Returns the norm of the projected step s
  //    s     -- The projected step upon return
  //    w     -- The direction vector w (unchanged)
  //    x     -- The anchor vector x (unchanged)
  //    alpha -- The step size (unchanged)
  //    model -- Contains the bound constraint information
  Real dgpstep(Vector<Real> &s, const Vector<Real> &w,
         const Vector<Real> &x, const Real alpha,
               TrustRegionModel<Real> &model) const {
    s.set(x); s.axpy(alpha,w);
    model.getBoundConstraint()->project(s);
    s.axpy(static_cast<Real>(-1),x);
    return s.norm();
  }

  // Compute minimal and maximal break points of x+alpha*s
  // with in the interval [xl,xu] specified by the bound constraint
  //   x     -- The anchor vector x (unchanged)
  //   s     -- The descent vector s (unchanged)
  //   model -- Contains the bound constraint information
  //   bpmin -- The minimum break point
  //   bpmax -- The maximum break point
  //   pwa   -- A primal working vector
  void dbreakpt(const Vector<Real> &x, const Vector<Real> &s,
                TrustRegionModel<Real> &model,
                Real &bpmin, Real &bpmax,
                Vector<Real> &pwa) {
    const Real zero(0), one(1);
    bpmin = one; bpmax = zero;
    Real lbpmin = one, lbpmax = zero, ubpmin = one, ubpmax = zero; 
    // Compute lower break points
    if (model.getBoundConstraint()->isLowerActivated()) {
      pwa.set(x);
      pwa.axpy(-one,*model.getBoundConstraint()->getLowerBound());
      pwa.applyBinary(lbp_,s);
      if (pwa.norm() != zero) {
        lbpmin = pwa.reduce(pmin_); 
        lbpmax = pwa.reduce(pmax_);
      }
    }
    // Compute upper break points
    if (model.getBoundConstraint()->isUpperActivated()) {
      pwa.set(*model.getBoundConstraint()->getUpperBound());
      pwa.axpy(-one,x);
      pwa.applyBinary(ubp_,s);
      if (pwa.norm() != zero) {
        ubpmin = pwa.reduce(pmin_); 
        ubpmax = pwa.reduce(pmax_);
      }
    }
    bpmin = std::min(lbpmin,ubpmin);
    bpmax = std::max(lbpmax,ubpmax);
    if (bpmin > bpmax) {
      bpmin = zero;
      bpmax = zero;
    }
    if (verbosity_ > 0) {
      std::cout << std::endl;
      std::cout << "  Computation of break points"          << std::endl;
      std::cout << "    Minimum break point:              " << bpmin      << std::endl;
      std::cout << "    Maximum break point:              " << bpmax      << std::endl;
    }
  }

  // Compute Cauchy point, i.e., the minimizer of q(P(x - alpha*g)-x)
  // subject to the trust region constraint ||P(x - alpha*g)-x|| <= del
  //   s     -- The Cauchy point upon return
  //   alpha -- The step length for the Cauchy point upon return
  //   x     -- The anchor vector x (unchanged)
  //   g     -- The (dual) gradient vector g (unchanged)
  //   del   -- The trust region radius (unchanged)
  //   model -- Contains the objective and bound constraint information
  //   pwa1  -- Primal working array
  //   pwa2  -- Primal working array
  //   dwa   -- Dual working array
  Real dcauchy(Vector<Real> &s, Real &alpha,
               const Vector<Real> &x, const Vector<Real> &g,
               const Real del, TrustRegionModel<Real> &model,
               Vector<Real> &pwa1, Vector<Real> &pwa2, Vector<Real> &dwa) {
    const Real half(0.5), one(1), mu0(0.01), interpf(0.1), extrapf(10);
    // const Real zero(0); // Unused
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    bool interp = false;
    Real q(0), gs(0), bpmin(0), bpmax(0), snorm(0);
    // Compute minimal and maximal break points of x[0] - alpha g[0]
    pwa1.set(g); pwa1.scale(-one);
    dbreakpt(x,pwa1,model,bpmin,bpmax,pwa2);
    // Compute s = P(x[0] - alpha g[0].dual())
    snorm = dgpstep(s,g,x,-alpha,model);
    if (snorm > del) {
      interp = true;
    }
    else {
      dynamic_cast<LinMoreModel<Real>&>(model).applyFullHessian(dwa,s,tol);
      gs = s.dot(g);
      q  = half * s.dot(dwa.dual()) + gs;
      interp = (q > mu0*gs);
    }
    // Either increase or decrease alpha to find approximate Cauchy point
    if (interp) {
      bool search = true;
      while (search) {
        alpha *= interpf;
        snorm = dgpstep(s,g,x,-alpha,model);
        if (snorm <= del) {
          dynamic_cast<LinMoreModel<Real>&>(model).applyFullHessian(dwa,s,tol);
          gs = s.dot(g);
          q  = half * s.dot(dwa.dual()) + gs;
          search = (q > mu0*gs);
        }
      }
    }
    else {
      bool search = true;
      Real alphas = alpha;
      while (search && alpha <= bpmax) {
        alpha *= extrapf;
        snorm = dgpstep(s,g,x,-alpha,model);
        if (snorm <= del) {
          dynamic_cast<LinMoreModel<Real>&>(model).applyFullHessian(dwa,s,tol);
          gs = s.dot(g);
          q  = half * s.dot(dwa.dual()) + gs;
          if (q < mu0*gs) {
            search = true;
            alphas = alpha;
          }
        }
        else {
          search = false;
        }
      }
      alpha = alphas;
      snorm = dgpstep(s,g,x,-alpha,model);
    }
    return snorm;
  }

  // Perform projected search to determine beta such that
  // q(P(x + beta*s)-x) <= mu0*g'(P(x + beta*s)-x) for mu0 in (0,1)
  //   x     -- The anchor vector x, upon return x = P(x + beta*s)
  //   s     -- The direction vector s, upon return s = P(x + beta*s) - x
  //   g     -- The free components of the gradient vector g (unchanged)
  //   model -- Contains objective and bound constraint information
  //   pwa   -- Primal working array
  //   dwa   -- Dual working array
  Real dprsrch(Vector<Real> &x, Vector<Real> &s,
               const Vector<Real> &g, TrustRegionModel<Real> &model,
               Vector<Real> &pwa, Vector<Real> &dwa) {
    const Real half(0.5), one(1), mu0(0.01), interpf(0.5);
    Real tol = std::sqrt(ROL_EPSILON<Real>());
    Real beta(1), snorm(0), q(0), gs(0), bpmin(0), bpmax(0);
    int nsteps = 0;
    // Compute break points of x+beta*s;
    dbreakpt(x,s,model,bpmin,bpmax,pwa);
    // Reduce beta until sufficient decrease is satisfied
    bool search = true;
    while (search && beta > bpmin) {
      nsteps++;
      snorm = dgpstep(pwa,s,x,beta,model);
      dynamic_cast<LinMoreModel<Real>&>(model).applyFreeHessian(dwa,pwa,x,tol);
      gs = pwa.dot(g);
      q  = half * s.dot(dwa.dual()) + gs;
      if (q <= mu0*gs) {
        search = false;
      }
      else {
        beta *= interpf;
      }
    }
    if (beta < one && beta < bpmin) {
      beta = bpmin;
    }
    snorm = dgpstep(pwa,s,x,beta,model);
    s.set(pwa);
    x.plus(s);
    if (verbosity_ > 0) {
      std::cout << std::endl;
      std::cout << "  Projected search"                     << std::endl;
      std::cout << "    Step length (beta):               " << beta       << std::endl;
    }
    return snorm;
  }

  // Compute sigma such that ||x+sigma*p||_inv(M) = del.  This is called
  // if dtrpcg detects negative curvature or if the step violates
  // the trust region bound
  //   xtx -- The dot product <x, inv(M)x> (unchanged)
  //   ptp -- The dot product <p, inv(M)p> (unchanged)
  //   ptx -- The dot product <p, inv(M)x> (unchanged)
  //   del -- Trust region radius (unchanged)
  Real dtrqsol(const Real xtx, const Real ptp, const Real ptx, const Real del) {
    const Real zero(0);
    Real dsq = del*del;
    Real rad = ptx*ptx + ptp*(dsq-xtx);
    rad = std::sqrt(std::max(rad,zero));
    Real sigma(0);
    if (ptx > zero) {
      sigma = (dsq-xtx)/(ptx+rad);
    }
    else if (rad > zero) {
      sigma = (rad-ptx)/ptp;
    }
    else {
      sigma = zero;
    }
    return sigma;
  }

  // Solve the trust region subproblem: minimize q(w) subject to the
  // trust region constraint ||w||_inv(M) <= del using the Steihaug-Toint
  // Conjugate Gradients algorithm
  //   w       -- The step vector to be computed
  //   iflag   -- Termination flag
  //   iter    -- Number of CG iterations
  //   del     -- Trust region radius (unchanged)
  //   model   -- Contains the objective and bound constraint information
  //   tol     -- Residual stopping tolerance (unchanged)
  //   stol    -- Preconditioned residual stopping tolerance (unchanged)
  //   itermax -- Maximum number of iterations
  //   p       -- Primal working array that stores the CG step
  //   q       -- Dual working array that stores the Hessian applied to p
  //   r       -- Primal working array that stores the preconditioned residual
  //   t       -- Dual working array that stores the residual
  Real dtrpcg(Vector<Real> &w, int &iflag, int &iter,
              const Vector<Real> &g, const Vector<Real> &x,
              const Real del, TrustRegionModel<Real> &model,
              const Real tol, const Real stol, const int itermax,
              Vector<Real> &p, Vector<Real> &q, Vector<Real> &r,
              Vector<Real> &t) {
    Real tol0 = std::sqrt(ROL_EPSILON<Real>());
    const Real zero(0), one(1), two(2);
    Real rho(0), tnorm(0), rnorm(0), rnorm0(0), kappa(0), beta(0), sigma(0), alpha(0), rtr(0);
    Real sMs(0), pMp(0), sMp(0);
    iter = 0; iflag = 0;
    // Initialize step
    w.zero();
    // Compute residual
    t.set(g); t.scale(-one);
    // Preconditioned residual
    dynamic_cast<LinMoreModel<Real>&>(model).applyFreePrecond(r,t,x,tol0);
    rho    = r.dot(t.dual());
    rnorm0 = std::sqrt(rho);
    if ( rnorm0 == zero ) {
      return zero;
    }
    // Initialize direction
    p.set(r);
    pMp = rho;
    // Iterate CG
    for (iter = 0; iter < itermax; ++iter) {
      // Apply Hessian to direction dir
      dynamic_cast<LinMoreModel<Real>&>(model).applyFreeHessian(q,p,x,tol0);
      // Compute sigma such that ||s+sigma*dir|| = del
      kappa = p.dot(q.dual());
      alpha = (kappa>zero) ? rho/kappa : zero;
      sigma = dtrqsol(sMs,pMp,sMp,del);
      // Check for negative curvature or if iterate exceeds trust region
      if (kappa <= zero || alpha >= sigma) {
        w.axpy(sigma,p);
        iflag = (kappa<=zero) ? 2 : 3;
        break;
      }
      // Update iterate and residuals
      w.axpy(alpha,p);
      t.axpy(-alpha,q);
      dynamic_cast<LinMoreModel<Real>&>(model).applyFreePrecond(r,t,x,tol0);
      // Exit if residual tolerance is met
      rtr   = r.dot(t.dual());
      rnorm = std::sqrt(rtr);
      tnorm = t.norm();
      if (rnorm <= stol || tnorm <= tol) {
        iflag = 0;
        break;
      }
      // Compute p = r + beta * p
      beta = rtr/rho;
      p.scale(beta); p.plus(r);
      rho  = rtr;
      // Update dot products
      //   sMs = <s, inv(M)s>
      //   sMp = <s, inv(M)p>
      //   pMp = <p, inv(M)p>
      sMs  = sMs + two*alpha*sMp + alpha*alpha*pMp;
      sMp  = beta*(sMp + alpha*pMp);
      pMp  = rho + beta*beta*pMp;
    }
    // Check iteration count
    if (iter == itermax) {
      iflag = 1;
    }
    if (iflag != 1) { 
      iter++;
    }
    return w.norm();
  }

};

}

#endif
