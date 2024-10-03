// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_FDIVERGENCE_HPP
#define ROL_FDIVERGENCE_HPP

#include "ROL_RandVarFunctional.hpp"
#include "ROL_Types.hpp"

/** @ingroup risk_group
    \class ROL::FDivergence
    \brief Provides a general interface for the F-divergence distributionally robust
    expectation.

    This class defines a risk measure \f$\mathcal{R}\f$ which arises in distributionally
    robust stochastic programming.  \f$\mathcal{R}\f$ is given by
    \f[
       \mathcal{R}(X) = \sup_{\vartheta\in\mathfrak{A}}
           \mathbb{E}[\vartheta X]
    \f]
    where \f$\mathfrak{A}\f$ is called the ambiguity (or uncertainty) set and 
    is defined by a constraint on the F-divergence, i.e.,
    \f[
       \mathfrak{A} = \{\vartheta\in\mathcal{X}^*\,:\,
         \mathbb{E}[\vartheta] = 1,\; \vartheta \ge 0,\;\text{and}\;
         \mathbb{E}[F(\vartheta)] \le \epsilon\}
    \f]
    where \f$F:\mathbb{R}\to[0,\infty]\f$ convex, lower semicontinuous and satisfies
    \f$F(1) = 1\f$ and \f$F(x) = \infty\f$ for \f$x < 0\f$.
    \f$\mathcal{R}\f$ is a law-invariant, coherent risk measure.  Moreover, by a
    duality argument, \f$\mathcal{R}\f$ can be reformulated as
    \f[
       \mathcal{R}(X) = \inf_{\lambda > 0,\,\mu}\left\{
             \lambda \epsilon + \mu + \mathbb{E}\left[
                (\lambda F)^*(X-\mu)\right]\right\}.
    \f]
    Here, \f$(\lambda F)^*\f$ denotes the Legendre-Fenchel transformation of
    \f$(\lambda F)\f$.
    ROL implements this by augmenting the optimization vector \f$x_0\f$ with
    the parameter \f$(\lambda,\mu)\f$, then minimizes jointly for
    \f$(x_0,\lambda,\mu)\f$.
*/

namespace ROL {

template<class Real>
class FDivergence : public RandVarFunctional<Real> {
private:
  Real thresh_;

  Real valLam_;
  Real valLam2_;
  Real valMu_;
  Real valMu2_;

  using RandVarFunctional<Real>::val_;
  using RandVarFunctional<Real>::gv_;
  using RandVarFunctional<Real>::g_;
  using RandVarFunctional<Real>::hv_;
  using RandVarFunctional<Real>::dualVector_;

  using RandVarFunctional<Real>::point_;
  using RandVarFunctional<Real>::weight_;

  using RandVarFunctional<Real>::computeValue;
  using RandVarFunctional<Real>::computeGradient;
  using RandVarFunctional<Real>::computeGradVec;
  using RandVarFunctional<Real>::computeHessVec;

  void checkInputs(void) const {
    Real zero(0);
    ROL_TEST_FOR_EXCEPTION((thresh_ <= zero), std::invalid_argument,
      ">>> ERROR (ROL::FDivergence): Threshold must be positive!");
  }

public:
  /** \brief Constructor.

      @param[in]     eps    is the tolerance for the F-divergence constraint
  */
  FDivergence(const Real thresh) : RandVarFunctional<Real>(), thresh_(thresh),
    valLam_(0),valLam2_(0), valMu_(0), valMu2_(0) {
    checkInputs();
  }

  /** \brief Constructor.

      @param[in]     parlist is a parameter list specifying inputs

      parlist should contain sublists "SOL"->"Risk Measure"->"F-Divergence" and
      within the "F-Divergence" sublist should have the following parameters
      \li "Threshold" (greater than 0)
  */
  FDivergence(ROL::ParameterList &parlist) : RandVarFunctional<Real>(),
    valLam_(0),valLam2_(0), valMu_(0), valMu2_(0) {
    ROL::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("F-Divergence");
    thresh_ = list.get<Real>("Threshold");
    checkInputs();
  }

  /** \brief Implementation of the scalar primal F function.

      @param[in]     x     is a scalar input
      @param[in]     deriv is the derivative order

      Upon return, Fprimal returns \f$F(x)\f$ or a derivative of \f$F(x)\f$.
  */
  virtual Real Fprimal(Real x, int deriv = 0) const = 0;

  /** \brief Implementation of the scalar dual F function.

      @param[in]     x     is a scalar input
      @param[in]     deriv is the derivative order

      Upon return, Fdual returns \f$F^*(x)\f$ or a derivative of \f$F^*(x)\f$.
      Here, \f$F^*\f$ denotes the Legendre-Fenchel transformation of \f$F\f$,
      i.e.,
      \f[
          F^*(y) = \sup_{x\in\mathbb{R}}\{xy - F(x)\}.
      \f]
  */
  virtual Real Fdual(Real x, int deriv = 0) const = 0;

  bool check(std::ostream &outStream = std::cout) const {
    const Real tol(std::sqrt(ROL_EPSILON<Real>()));
    bool flag = true;

    Real x  = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
    Real t  = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
    Real fp = Fprimal(x);
    Real fd = Fdual(t);
    outStream << "Check Fenchel-Young Inequality: F(x) + F*(t) >= xt" << std::endl;
    outStream << "x       = " << x                << std::endl;
    outStream << "t       = " << t                << std::endl;
    outStream << "F(x)    = " << fp               << std::endl;
    outStream << "F*(t)   = " << fd               << std::endl;
    outStream << "Is Valid? " << (fp+fd >= x*t)   << std::endl;
    flag = (fp+fd >= x*t) ? flag : false;

    x  = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
    t  = Fprimal(x,1);
    fp = Fprimal(x);
    fd = Fdual(t);
    outStream << "Check Fenchel-Young Equality: F(x) + F(t) = xt for t = d/dx F(x)" << std::endl;
    outStream << "x       = " << x                << std::endl;
    outStream << "t       = " << t                << std::endl;
    outStream << "F(x)    = " << fp               << std::endl;
    outStream << "F*(t)   = " << fd               << std::endl;
    outStream << "Is Valid? " << (std::abs(fp+fd - x*t)<=tol) << std::endl; 
    flag = (std::abs(fp+fd - x*t)<=tol) ? flag : false;

    t  = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
    x  = Fdual(t,1);
    fp = Fprimal(x);
    fd = Fdual(t);
    outStream << "Check Fenchel-Young Equality: F(x) + F(t) = xt for x = d/dt F*(t)" << std::endl;
    outStream << "x       = " << x                << std::endl;
    outStream << "t       = " << t                << std::endl;
    outStream << "F(x)    = " << fp               << std::endl;
    outStream << "F*(t)   = " << fd               << std::endl;
    outStream << "Is Valid? " << (std::abs(fp+fd - x*t)<=tol) << std::endl; 
    flag = (std::abs(fp+fd - x*t)<=tol) ? flag : false;

    return flag;
  }

  void initialize(const Vector<Real> &x) {
    RandVarFunctional<Real>::initialize(x);
    valLam_ = 0; valLam2_ = 0; valMu_ = 0; valMu2_ = 0;
  }

  // Value update and get functions
  void updateValue(Objective<Real>         &obj,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   Real                    &tol) {
    Real val  = computeValue(obj,x,tol);
    Real xlam = xstat[0];
    Real xmu  = xstat[1];
    Real r    = Fdual((val-xmu)/xlam,0);
    val_     += weight_ * r;
  }

  Real getValue(const Vector<Real>      &x,
                const std::vector<Real> &xstat,
                SampleGenerator<Real>   &sampler) {
    Real val(0);
    sampler.sumAll(&val_,&val,1);
    Real xlam = xstat[0];
    Real xmu  = xstat[1];
    return xlam*(thresh_ + val) + xmu;
  }

  // Gradient update and get functions
  void updateGradient(Objective<Real>         &obj,
                      const Vector<Real>      &x,
                      const std::vector<Real> &xstat,
                      Real                    &tol) {
    Real val  = computeValue(obj,x,tol);
    Real xlam = xstat[0];
    Real xmu  = xstat[1];
    Real inp  = (val-xmu)/xlam;
    Real r0 = Fdual(inp,0), r1 = Fdual(inp,1);

    if (std::abs(r0) >= ROL_EPSILON<Real>()) {
      val_ += weight_ * r0;
    }
    if (std::abs(r1) >= ROL_EPSILON<Real>()) {
      valLam_ -= weight_ * r1 * inp;
      valMu_  -= weight_ * r1;
      computeGradient(*dualVector_,obj,x,tol);
      g_->axpy(weight_*r1,*dualVector_);
    }
  }

  void getGradient(Vector<Real>            &g,
                   std::vector<Real>       &gstat,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   SampleGenerator<Real>   &sampler) {
    std::vector<Real> mygval(3), gval(3);
    mygval[0] = val_;
    mygval[1] = valLam_;
    mygval[2] = valMu_;
    sampler.sumAll(&mygval[0],&gval[0],3);

    gstat[0] = thresh_ + gval[0] + gval[1];
    gstat[1] = static_cast<Real>(1) + gval[2];

    sampler.sumAll(*g_,g);
  }

  void updateHessVec(Objective<Real>         &obj,
                     const Vector<Real>      &v,
                     const std::vector<Real> &vstat,
                     const Vector<Real>      &x,
                     const std::vector<Real> &xstat,
                     Real                    &tol) {
    Real val  = computeValue(obj,x,tol);
    Real xlam = xstat[0];
    Real xmu  = xstat[1];
    Real vlam = vstat[0];
    Real vmu  = vstat[1];
    Real inp  = (val-xmu)/xlam;
    Real r1 = Fdual(inp,1), r2 = Fdual(inp,2);
    if (std::abs(r2) >= ROL_EPSILON<Real>()) {
      Real gv   = computeGradVec(*dualVector_,obj,v,x,tol);
      val_     += weight_ * r2 * inp;
      valLam_  += weight_ * r2 * inp * inp;
      valLam2_ -= weight_ * r2 * gv * inp;
      valMu_   += weight_ * r2;
      valMu2_  -= weight_ * r2 * gv;
      hv_->axpy(weight_ * r2 * (gv - vmu - vlam*inp)/xlam, *dualVector_);
    }
    if (std::abs(r1) >= ROL_EPSILON<Real>()) {
      computeHessVec(*dualVector_,obj,v,x,tol);
      hv_->axpy(weight_ * r1, *dualVector_);
    }
  }

  void getHessVec(Vector<Real>            &hv,
                  std::vector<Real>       &hvstat,
                  const Vector<Real>      &v,
                  const std::vector<Real> &vstat,
                  const Vector<Real>      &x,
                  const std::vector<Real> &xstat,
                  SampleGenerator<Real>   &sampler) {
    std::vector<Real> myhval(5), hval(5);
    myhval[0] = val_;
    myhval[1] = valLam_;
    myhval[2] = valLam2_;
    myhval[3] = valMu_;
    myhval[4] = valMu2_;
    sampler.sumAll(&myhval[0],&hval[0],5);

    std::vector<Real> stat(2);
    Real xlam = xstat[0];
    //Real xmu  = xstat[1];
    Real vlam = vstat[0];
    Real vmu  = vstat[1];
    hvstat[0] = (vlam * hval[1] + vmu * hval[0] + hval[2])/xlam;
    hvstat[1] = (vlam * hval[0] + vmu * hval[3] + hval[4])/xlam;

    sampler.sumAll(*hv_,hv);
  }
};

}

#endif
