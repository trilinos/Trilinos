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

#ifndef ROL_REDUCEDDYNAMICOBJECTIVE_HPP
#define ROL_REDUCEDDYNAMICOBJECTIVE_HPP

#include "ROL_Ptr.hpp"
#include "ROL_Sketch.hpp"
#include "ROL_Objective.hpp"
#include "ROL_DynamicObjective.hpp"
#include "ROL_DynamicConstraint.hpp"


/** @ingroup func_group
    \class ROL::ReducedDynamicObjective
    \brief Defines the reduced time-dependent objective function interface
           for simulation-based optimization.

    This objective function implements the implicitly defined objective
    function given by
    \f[
       F(z) := g\left(\sum_{n=1}^{N_t} f_n(u_{n-1},u_n(z),z_n)\right)
    \f]
    where \f$g:\mathbb{R}\to\mathbb{R}\f$,
    \f$f_n:\mathcal{U}\times\mathcal{U}\times\mathcal{Z}\to\mathbb{R}\f$,
    and \f$u_n\in\mathcal{U}\f$ solves the system of equations
    \f[
       c_n(u_{n-1},u_n,z_n) = 0,\quad n=1,\ldots,N_t
    \f]
    with \f$u_0\f$ provided.
*/


namespace ROL {

template<typename Real> 
class ReducedDynamicObjective : public Objective<Real> {
  using size_type = typename std::vector<Real>::size_type;
private:
  const Ptr<DynamicObjective<Real>>  obj_;
  const Ptr<DynamicConstraint<Real>> con_;
  const Ptr<Vector<Real>>            u0_;
  const std::vector<TimeStamp<Real>> timeStamp_;
  const size_type                    Nt_;
  const bool                         useSketch_;
  const int                          rank_;
  Ptr<Sketch<Real>>                  sketch_;
  std::vector<Ptr<Vector<Real>>>     uhist_;
  std::vector<Ptr<Vector<Real>>>     lhist_;
  Ptr<Vector<Real>>                  cprimal_;
  Ptr<Vector<Real>>                  udual_;
  Ptr<Vector<Real>>                  rhs_;
  Ptr<Vector<Real>>                  zdual_;
  Real                               val_;
  bool                               isValueComputed_;

  PartitionedVector<Real>& partition ( Vector<Real>& x ) const {
    return static_cast<PartitionedVector<Real>&>(x);
  }

  const PartitionedVector<Real>& partition ( const Vector<Real>& x ) const {
    return static_cast<const PartitionedVector<Real>&>(x);
  }

  void setTerminalCondition(Vector<Real> &l,
                      const Vector<Real> &uold, const Vector<Real> &unew,
                      const Vector<Real> &z,    const TimeStamp<Real> &ts) {
    obj_->gradient_un(*udual_, uold, unew, z, ts);
    con_->applyInverseAdjointJacobian_un(l, *udual_, uold, unew, z, ts);
  }

  void updateGradient(Vector<Real> &g,    const Vector<Real> &l,
                const Vector<Real> &uold, const Vector<Real> &unew,
                const Vector<Real> &z,    const TimeStamp<Real> &ts) {
    const Real one(1);
    obj_->gradient_z(g, uold, unew, z, ts);
    con_->applyAdjointJacobian_z(*zdual_, l, uold, unew, z, ts);
    g.axpy(-one,*zdual_);
  }

  void computeAdjointRHS(Vector<Real> &rhs,  const Vector<Real>    &l,
                   const Vector<Real> &uold, const Vector<Real>    &unew,
                   const Vector<Real> &z,    const TimeStamp<Real> &ts) {
    const Real one(1);
    obj_->gradient_uo(rhs, uold, unew, z, ts);
    con_->applyAdjointJacobian_uo(*udual_, l, uold, unew, z, ts);
    rhs.axpy(-one,*udual_);
  }

  void advanceAdjoint(Vector<Real> &l,          Vector<Real>    &rhs,
                const Vector<Real> &uold, const Vector<Real>    &unew,
                const Vector<Real> &z,    const TimeStamp<Real> &ts) {
    obj_->gradient_un(*udual_, uold, unew, z, ts);
    rhs.plus(*udual_);
    con_->applyInverseAdjointJacobian_un(l, rhs, uold, unew, z, ts);
  }

public:
  ReducedDynamicObjective(const Ptr<DynamicObjective<Real>>  &obj,
                          const Ptr<DynamicConstraint<Real>> &con,
                          const Ptr<Vector<Real>>            &u0,
                          const Ptr<Vector<Real>>            &zvec,
                          const Ptr<Vector<Real>>            &cvec,
                          const std::vector<TimeStamp<Real>> &timeStamp,
                          const bool                          useSketch = false,
                          const int                           rank      = 10)
    : obj_(obj), con_(con), u0_(u0), timeStamp_(timeStamp),
      Nt_(timeStamp.size()), useSketch_(useSketch), rank_(rank),
      sketch_(nullPtr), isValueComputed_(false) {
    uhist_.clear(); lhist_.clear();
    if (useSketch_) { // Only maintain a sketch of the state time history
      sketch_ = makePtr<Sketch<Real>>(*u0_,static_cast<int>(Nt_),rank_);
      uhist_.push_back(u0_->clone());
      uhist_.push_back(u0_->clone());
      lhist_.push_back(cvec->dual().clone());
    }
    else {            // Store entire state time history
      for (size_type i = 0; i < Nt_; ++i) {
        uhist_.push_back(u0_->clone());
        lhist_.push_back(cvec->dual().clone());
      }
    }
    cprimal_ = cvec->clone();
    udual_   = u0_->dual().clone();
    rhs_     = u0_->dual().clone();
    zdual_   = zvec->dual().clone();
  }

  virtual Real outerFunction(const Real x, const int deriv = 0) const {
    if (deriv < 0 || deriv > 2) {
      throw Exception::NotImplemented(">>> ROL::ReducedDynamicConstraint::outFunction: deriv must be 0, 1, or 2!");
    }
    const Real zero(0), one(1);
    return (deriv == 0 ? x : (deriv == 1 ? one : zero));
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    if (useSketch_) {
      sketch_->update();
    }
    for (size_type i = 0; i < uhist_.size(); ++i) {
      uhist_[i]->zero();
    }
    val_ = static_cast<Real>(0);
    isValueComputed_ = false;
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    if (!isValueComputed_) {
      const Real one(1);
      const PartitionedVector<Real> &xp = partition(x);
      // Set initial condition
      uhist_[0]->set(*u0_);
      if (useSketch_) {
        sketch_->advance(one,*uhist_[0],0,one);
      }
      // Run time stepper
      Real valk(0);
      size_type index;
      for (size_type k = 1; k < Nt_; ++k) {
        index = (useSketch_ ? 1 : k);
        // Solve state on current time interval
        con_->solve(*cprimal_, *uhist_[index-1], *uhist_[index], *xp.get(k), timeStamp_[k]);
        // Compute objective function value on current time interval
        valk  = obj_->value(*uhist_[index-1], *uhist_[index], *xp.get(k), timeStamp_[k]);
        // Update total objective function value
        val_ += valk;
        // Sketch state
        if (useSketch_) {
          sketch_->advance(one,*uhist_[index],static_cast<int>(k),one);
          uhist_[index-1]->set(*uhist_[index]);
        }
      }
      isValueComputed_ = true;
    }
    return outerFunction(val_,0);
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    PartitionedVector<Real>       &gp = partition(g);
    const PartitionedVector<Real> &xp = partition(x);
    const Real one(1);
    size_type uindex = (useSketch_ ? 1 : Nt_-1);
    size_type lindex = (useSketch_ ? 0 : Nt_-1);
    // Must first compute the value
    if (!isValueComputed_) {
      Real val = value(x,tol);
    }
    Real gval = outerFunction(val_,1);
    // Recover state from sketch
    if (useSketch_) {
      uhist_[1]->set(*uhist_[0]);
      sketch_->reconstruct(*uhist_[0],Nt_-2);
    }
    // Solve for terminal condition
    setTerminalCondition(*lhist_[lindex],
                         *uhist_[uindex-1], *uhist_[uindex],
                         *xp.get(Nt_-1),    timeStamp_[Nt_-1]);
    // Update gradient on terminal interval
    updateGradient(*gp.get(Nt_-1),    *lhist_[lindex],
                   *uhist_[uindex-1], *uhist_[uindex],
                   *xp.get(Nt_-1),    timeStamp_[Nt_-1]);
    // Run reverse time stepper
    for (size_type k = Nt_-2; k > 0; --k) {
      // Compute k+1 component of rhs
      computeAdjointRHS(*rhs_,             *lhist_[lindex],
                        *uhist_[uindex-1], *uhist_[uindex],
                        *xp.get(k+1),      timeStamp_[k+1]);
      // Recover state from sketch
      if (useSketch_) {
        uhist_[1]->set(*uhist_[0]);
        if (k==1) {
          uhist_[0]->set(*u0_);
        }
        else {
          sketch_->reconstruct(*uhist_[0],k-1);
        }
      }
      uindex = (useSketch_ ? 1 : k);
      lindex = (useSketch_ ? 0 : k);
      // Solve for adjoint on interval k
      advanceAdjoint(*lhist_[lindex],   *rhs_,
                     *uhist_[uindex-1], *uhist_[uindex],
                     *xp.get(k),        timeStamp_[k]);
      // Update gradient on interval k
      updateGradient(*gp.get(k),        *lhist_[lindex],
                     *uhist_[uindex-1], *uhist_[uindex],
                     *xp.get(k),        timeStamp_[k]);
    }
    g.scale(gval);
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    Objective<Real>::hessVec(hv,v,x,tol);
  }
};

} // namespace ROL


#endif // ROL_REDUCEDDYNAMICOBJECTIVE_HPP

