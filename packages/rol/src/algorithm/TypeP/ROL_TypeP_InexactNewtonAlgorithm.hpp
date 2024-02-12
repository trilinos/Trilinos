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

#ifndef ROL_TYPEP_INEXACTNEWTONALGORITHM_HPP
#define ROL_TYPEP_INEXACTNEWTONALGORITHM_HPP

#include "ROL_TypeP_Algorithm.hpp"

/** \class ROL::TypeP::InexactNewtonAlgorithm
    \brief Provides an interface to run the inexact proximal Newton algorithm.
*/

namespace ROL {
namespace TypeP {

template<typename Real>
class InexactNewtonAlgorithm : public TypeP::Algorithm<Real> {
private:
  int t0_;
  bool initProx_;

  int maxit_;         ///< Maximum number of line search steps (default: 20)
  Real rhodec_;       ///< Backtracking rate (default: 0.5)
  Real c1_;           ///< Sufficient Decrease Parameter (default: 1e-4)
  Real sigma1_;       ///< Lower safeguard for quadratic line search (default: 0.1)
  Real sigma2_;       ///< Upper safeguard for quadratic line search (default: 0.9)
  Real sp_tol1_;
  Real sp_tol2_;
  Real sp_exp_;
  Real sp_tol_min_;
  std::string algoName_;

  ParameterList list_;

  int ls_nfval_, spgIter_, nhess_;
  int verbosity_;
  bool writeHeader_;

  using TypeP::Algorithm<Real>::pgstep;
  using TypeP::Algorithm<Real>::status_;
  using TypeP::Algorithm<Real>::state_;

  class NewtonObj : public Objective<Real> {
  private:
    const Ptr<Objective<Real>> obj_;
    const Ptr<Vector<Real>>    x_;
    const Ptr<Vector<Real>>    g_;
    Ptr<Vector<Real>>          pwa_, dwa_, Hx_;
    bool                       isHessApplied_;
    int                        nhess_;

  public:
    NewtonObj(const Ptr<Objective<Real>> &obj, const Vector<Real> &x, const Vector<Real> &g)
      : obj_(obj), x_(x.clone()), g_(g.clone()), pwa_(x.clone()),
        dwa_(g.clone()), Hx_(g.clone()), isHessApplied_(false), nhess_(0) {}
    void update(const Vector<Real> &x, UpdateType type, int iter) {
      isHessApplied_ = false;
    }
    Real value(const Vector<Real> &x, Real &tol) {
      const Real half(0.5), one(1);
      pwa_->set(x); pwa_->axpy(-one,*x_);
      if (!isHessApplied_) {
        obj_->hessVec(*Hx_,*pwa_,*x_,tol); nhess_++;
        isHessApplied_ = true;
      }
      dwa_->set(*Hx_);
      dwa_->scale(half);
      dwa_->plus(*g_);
      return dwa_->apply(*pwa_);
    }
    void gradient(Vector<Real> &g, const Vector<Real> &x, Real &tol) {
      const Real one(1);
      if (!isHessApplied_) {
        pwa_->set(x); pwa_->axpy(-one,*x_);
        obj_->hessVec(*Hx_,*pwa_,*x_,tol); nhess_++;
        isHessApplied_ = true;
      }
      g.set(*Hx_);
      g.plus(*g_);
    }
    void hessVec(Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol) {
      obj_->hessVec(hv,v,*x_,tol); nhess_++;
    }
    int numHessVec(bool reset = true) {
      int nhess = nhess_;
      if (reset) nhess_ = 0;
      return nhess;
    }
    void setData(const Vector<Real> &x, const Vector<Real> &g) {
      x_->set(x);
      g_->set(g);
    }
  };

  void initialize(Vector<Real>          &x,
                  const Vector<Real>    &g,
                  Objective<Real>       &sobj,
                  Objective<Real>       &nobj,
                  Vector<Real>          &dg,
                  Vector<Real>          &px,
                  std::ostream &outStream = std::cout); 

public:

  InexactNewtonAlgorithm(ParameterList &list);

  using TypeP::Algorithm<Real>::run;

  void run( Vector<Real>          &x,
            const Vector<Real>    &g, 
            Objective<Real>       &sobj,
            Objective<Real>       &nobj,
            std::ostream          &outStream = std::cout) override;

  void writeHeader( std::ostream& os ) const override;

  void writeName( std::ostream& os ) const override;

  void writeOutput( std::ostream &os, bool write_header = false ) const override;

}; // class ROL::TypeP::InexactNewtonAlgorithm

} // namespace TypeP
} // namespace ROL

#include "ROL_TypeP_InexactNewtonAlgorithm_Def.hpp"

#endif
