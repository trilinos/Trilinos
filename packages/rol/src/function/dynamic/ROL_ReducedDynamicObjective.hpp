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
       F(z) := \sum_{n=1}^{N_t} f_n(u_{n-1}(z),u_n(z),z_n)
    \f]
    where \f$f_n:\mathcal{U}\times\mathcal{U}\times\mathcal{Z}\to\mathbb{R}\f$,
    and \f$u_n\in\mathcal{U}\f$ solves the system of equations
    \f[
       c_n(u_{n-1},u_n,z_n) = 0,\quad n=1,\ldots,N_t
    \f]
    with \f$u_0\f$ provided.

    Disclaimer: This is currently only set up for single step time integrators
    and piecewise-constant-in-time controls.
*/


namespace ROL {

template<typename Real> 
class ReducedDynamicObjective : public Objective<Real> {
  using size_type = typename std::vector<Real>::size_type;
private:
  // Problem data.
  const Ptr<DynamicObjective<Real>>  obj_;
  const Ptr<DynamicConstraint<Real>> con_;
  const Ptr<Vector<Real>>            u0_;
  const std::vector<TimeStamp<Real>> timeStamp_;
  const size_type                    Nt_;
  // General sketch information.
  const bool                         useSketch_;
  // State sketch information.
  int                                rankState_;
  Ptr<Sketch<Real>>                  stateSketch_;
  // Adjoint sketch information.
  const int                          rankAdjoint_;
  Ptr<Sketch<Real>>                  adjointSketch_;
  // State sensitivity sketch information.
  const int                          rankStateSens_;
  Ptr<Sketch<Real>>                  stateSensSketch_;
  // Inexactness information.
  const bool                         useInexact_;
  const Real                         updateFactor_;
  const int                          maxRank_;
  // Vector storage for intermediate computations.
  std::vector<Ptr<Vector<Real>>>     uhist_;              // State history.
  std::vector<Ptr<Vector<Real>>>     lhist_;              // Adjoint history.
  std::vector<Ptr<Vector<Real>>>     whist_;              // State sensitivity history.
  std::vector<Ptr<Vector<Real>>>     phist_;              // Adjoint sensitivity history.
  Ptr<Vector<Real>>                  cprimal_;            // Primal constraint vector.
  Ptr<Vector<Real>>                  crhs_;               // Primal constraint vector.
  Ptr<Vector<Real>>                  udual_;              // Dual state vector.
  Ptr<Vector<Real>>                  rhs_;                // Dual state vector.
  Ptr<Vector<Real>>                  zdual_;              // Dual control vector.
  Real                               val_;                // Objective function value.
  // Flags.
  bool                               isValueComputed_;    // Whether objective has been evaluated.
  bool                               isStateComputed_;    // Whether state has been solved.
  bool                               isAdjointComputed_;  // Whether adjoint has been solved.
  bool                               useHessian_;         // Whether to use Hessian or FD.
  bool                               useSymHess_;         // Whether to use symmetric Hessian approximation.

  PartitionedVector<Real>& partition ( Vector<Real>& x ) const {
    return static_cast<PartitionedVector<Real>&>(x);
  }

  const PartitionedVector<Real>& partition ( const Vector<Real>& x ) const {
    return static_cast<const PartitionedVector<Real>&>(x);
  }

public:
  ReducedDynamicObjective(const Ptr<DynamicObjective<Real>>  &obj,
                          const Ptr<DynamicConstraint<Real>> &con,
                          const Ptr<Vector<Real>>            &u0,
                          const Ptr<Vector<Real>>            &zvec,
                          const Ptr<Vector<Real>>            &cvec,
                          const std::vector<TimeStamp<Real>> &timeStamp,
                          ROL::ParameterList                 &pl)
    : obj_               ( obj ),                                           // Dynamic objective function.
      con_               ( con ),                                           // Dynamic constraint function.
      u0_                ( u0 ),                                            // Initial condition.
      timeStamp_         ( timeStamp ),                                     // Vector of time stamps.
      Nt_                ( timeStamp.size() ),                              // Number of time intervals.
      useSketch_         ( pl.get("Use Sketching", true) ),                 // Use state sketch if true.
      rankState_         ( pl.get("State Rank", 10) ),                      // Rank of state sketch.
      stateSketch_       ( nullPtr ),                                       // State sketch object.
      rankAdjoint_       ( pl.get("Adjoint Rank", 10) ),                    // Rank of adjoint sketch.
      adjointSketch_     ( nullPtr ),                                       // Adjoint sketch object.
      rankStateSens_     ( pl.get("State Sensitivity Rank", 10) ),          // Rank of state sensitivity sketch.
      stateSensSketch_   ( nullPtr ),                                       // State sensitivity sketch object.
      useInexact_        ( pl.get("Adaptive Rank", false) ),                // Update rank adaptively.
      updateFactor_      ( pl.get("Rank Update Factor", 2.0) ),             // Rank update factor.
      maxRank_           ( pl.get("Maximum Rank", Nt_) ),                   // Maximum rank.
      isValueComputed_   ( false ),                                         // Flag indicating whether value has been computed.
      isStateComputed_   ( false ),                                         // Flag indicating whether state has been computed.
      isAdjointComputed_ ( false ),                                         // Flag indicating whether adjoint has been computed.
      useHessian_        ( pl.get("Use Hessian", true) ),                   // Flag indicating whether to use the Hessian.
      useSymHess_        ( pl.get("Use Only Sketched Sensitivity", true) ) {// Flag indicating whether to use symmetric sketched Hessian.
    uhist_.clear(); lhist_.clear(); whist_.clear(); phist_.clear();
    if (useSketch_) { // Only maintain a sketch of the state time history
      stateSketch_ = makePtr<Sketch<Real>>(*u0_,static_cast<int>(Nt_)-1,rankState_);
      uhist_.push_back(u0_->clone());
      uhist_.push_back(u0_->clone());
      lhist_.push_back(cvec->dual().clone());
      adjointSketch_ = makePtr<Sketch<Real>>(*u0_,static_cast<int>(Nt_)-1,rankAdjoint_);
      if (useHessian_) {
        stateSensSketch_ = makePtr<Sketch<Real>>(*u0_,static_cast<int>(Nt_)-1,rankStateSens_);
        whist_.push_back(u0_->clone());
        whist_.push_back(u0_->clone());
        phist_.push_back(cvec->dual().clone());
      }
    }
    else {            // Store entire state time history
      for (size_type i = 0; i < Nt_; ++i) {
        uhist_.push_back(u0_->clone());
        lhist_.push_back(cvec->dual().clone());
        if (useHessian_) {
          whist_.push_back(u0_->clone());
          phist_.push_back(cvec->dual().clone());
        }
      }
    }
    cprimal_ = cvec->clone();
    crhs_    = cvec->clone();
    udual_   = u0_->dual().clone();
    rhs_     = u0_->dual().clone();
    zdual_   = zvec->dual().clone();
  }

  Ptr<Vector<Real>> makeDynamicVector(const Vector<Real> &x) const {
    return ROL::PartitionedVector<Real>::create(x, Nt_);
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    if (useSketch_) {
      stateSketch_->update();
      if (flag == true) {
        adjointSketch_->update();
        if (useHessian_) {
          stateSensSketch_->update();
        }
      }
    }
    for (size_type i = 0; i < uhist_.size(); ++i) {
      uhist_[i]->zero();
    }
    val_ = static_cast<Real>(0);
    isValueComputed_ = false;
    isStateComputed_ = false;
    if (flag == true) {
      isAdjointComputed_ = false;
    }
  }

  Real value( const Vector<Real> &x, Real &tol ) {
    if (!isValueComputed_) {
      const Real one(1);
      const PartitionedVector<Real> &xp = partition(x);
      // Set initial condition
      uhist_[0]->set(*u0_);
      if (useSketch_) {
        stateSketch_->update();
        //stateSketch_->advance(one,*uhist_[0],0,one);
      }
      // Run time stepper
      Real valk(0);
      size_type index;
      for (size_type k = 1; k < Nt_; ++k) {
        index = (useSketch_ ? 1 : k);
        // Update dynamic constraint and objective
        con_->update(*uhist_[index-1], *uhist_[index], *xp.get(k), timeStamp_[k]);
        obj_->update(*uhist_[index-1], *uhist_[index], *xp.get(k), timeStamp_[k]);
        // Solve state on current time interval
        con_->solve(*cprimal_, *uhist_[index-1], *uhist_[index], *xp.get(k), timeStamp_[k]);
        // Compute objective function value on current time interval
        valk  = obj_->value(*uhist_[index-1], *uhist_[index], *xp.get(k), timeStamp_[k]);
        // Update total objective function value
        val_ += valk;
        // Sketch state
        if (useSketch_) {
          stateSketch_->advance(one,*uhist_[1],static_cast<int>(k)-1,one);
          uhist_[0]->set(*uhist_[1]);
        }
      }
      isValueComputed_ = true;
      isStateComputed_ = true;
    }
    return val_;
  }

  void gradient( Vector<Real> &g, const Vector<Real> &x, Real &tol ) {
    PartitionedVector<Real>       &gp = partition(g);
    gp.get(0)->zero(); // zero for the nonexistant zeroth control interval
    const PartitionedVector<Real> &xp = partition(x);
    const Real one(1);
    size_type uindex = (useSketch_ ? 1 : Nt_-1);
    size_type lindex = (useSketch_ ? 0 : Nt_-1);
    // Must first compute the value
    solveState(x);
    if (useSketch_ && useInexact_) {
      tol = updateSketch(x,tol);
    }
    // Recover state from sketch
    if (useSketch_) {
      uhist_[1]->set(*uhist_[0]);
      stateSketch_->reconstruct(*uhist_[0],static_cast<int>(Nt_)-3);
      if (isAdjointComputed_) {
        adjointSketch_->reconstruct(*lhist_[0],static_cast<int>(Nt_)-2);
      }
    }
    // Update dynamic constraint and objective
    con_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(Nt_-1), timeStamp_[Nt_-1]);
    obj_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(Nt_-1), timeStamp_[Nt_-1]);
    // Solve for terminal condition
    if (!isAdjointComputed_) {
      setTerminalCondition(*lhist_[lindex],
                           *uhist_[uindex-1], *uhist_[uindex],
                           *xp.get(Nt_-1),    timeStamp_[Nt_-1]);
      if (useSketch_) {
        adjointSketch_->advance(one,*lhist_[0],static_cast<int>(Nt_)-2,one);
      }
    }
    // Update gradient on terminal interval
    updateGradient(*gp.get(Nt_-1),    *lhist_[lindex],
                   *uhist_[uindex-1], *uhist_[uindex],
                   *xp.get(Nt_-1),    timeStamp_[Nt_-1]);
    // Run reverse time stepper
    for (size_type k = Nt_-2; k > 0; --k) {
      if (!isAdjointComputed_) {
        // Compute k+1 component of rhs
        computeAdjointRHS(*rhs_,             *lhist_[lindex],
                          *uhist_[uindex-1], *uhist_[uindex],
                          *xp.get(k+1),      timeStamp_[k+1]);
      }
      uindex = (useSketch_ ? 1 : k);
      lindex = (useSketch_ ? 0 : k);
      // Recover state from sketch
      if (useSketch_) {
        uhist_[1]->set(*uhist_[0]);
        if (k==1) {
          uhist_[0]->set(*u0_);
        }
        else {
          stateSketch_->reconstruct(*uhist_[0],static_cast<int>(k)-2);
        }
        if (isAdjointComputed_) {
          adjointSketch_->reconstruct(*lhist_[0],static_cast<int>(k)-1);
        }
      }
      // Update dynamic constraint and objective
      con_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(k), timeStamp_[k]);
      obj_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(k), timeStamp_[k]);
      // Solve for adjoint on interval k
      if (!isAdjointComputed_) {
        advanceAdjoint(*lhist_[lindex],   *rhs_,
                       *uhist_[uindex-1], *uhist_[uindex],
                       *xp.get(k),        timeStamp_[k]);
        if (useSketch_) {
          adjointSketch_->advance(one,*lhist_[0],static_cast<int>(k)-1,one);
        }
      }
      // Update gradient on interval k
      updateGradient(*gp.get(k),        *lhist_[lindex],
                     *uhist_[uindex-1], *uhist_[uindex],
                     *xp.get(k),        timeStamp_[k]);
    }
    isAdjointComputed_ = true;
  }

  void hessVec( Vector<Real> &hv, const Vector<Real> &v, const Vector<Real> &x, Real &tol ) {
    if (useHessian_) {
      // Must first solve the state and adjoint equations
      solveState(x);
      solveAdjoint(x);
      // Now compute Hessian
      const Real one(1);
      const PartitionedVector<Real> &xp  = partition(x);
      const PartitionedVector<Real> &vp  = partition(v);
      PartitionedVector<Real>       &hvp = partition(hv);
      hvp.get(0)->zero(); // zero for the nonexistant zeroth control interval
      // Compute state sensitivity
      whist_[0]->zero();
      if (useSketch_) {
        stateSensSketch_->update();
        //stateSensSketch_->advance(one,*whist_[0],0,one);
        uhist_[0]->set(*u0_);
      }
      size_type uindex, lindex;
      for (size_type k = 1; k < Nt_; ++k) {
        uindex = (useSketch_ ? 1 : k);
        lindex = (useSketch_ ? 0 : k);
        // Reconstruct sketched state
        if (useSketch_) {
          stateSketch_->reconstruct(*uhist_[1],static_cast<int>(k)-1);
          adjointSketch_->reconstruct(*lhist_[0],static_cast<int>(k)-1);
        }
        // Update dynamic constraint and objective
        con_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(k), timeStamp_[k]);
        obj_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(k), timeStamp_[k]);
        // Advance state sensitivity on current time interval
        advanceStateSens(*whist_[uindex],
                         *vp.get(k),        *whist_[uindex-1],
                         *uhist_[uindex-1], *uhist_[uindex],
                         *xp.get(k),        timeStamp_[k]);
        // Set control only Hessian
        computeControlHessLag(*hvp.get(k),
                              *vp.get(k),        *lhist_[lindex],
                              *uhist_[uindex-1], *uhist_[uindex],
                              *xp.get(k),        timeStamp_[k]);
        if (!useSymHess_) {
          // Add mixed derivative Hessian
          addMixedHessLag(*hvp.get(k),       *lhist_[lindex],
                          *whist_[uindex-1], *whist_[uindex],
                          *uhist_[uindex-1], *uhist_[uindex],
                          *xp.get(k),        timeStamp_[k]);
        }
        // Sketch state sensitivity
        if (useSketch_) {
          stateSensSketch_->advance(one,*whist_[1],static_cast<int>(k)-1,one);
          whist_[0]->set(*whist_[1]);
          uhist_[0]->set(*uhist_[1]);
        }
      }

      // Compute adjoint sensitivity
      uindex = (useSketch_ ? 1 : Nt_-1);
      lindex = (useSketch_ ? 0 : Nt_-1);
      if (useSketch_) {
        // Recover terminal state
        uhist_[1]->set(*uhist_[0]);
        stateSketch_->reconstruct(*uhist_[0],static_cast<int>(Nt_)-3);
        // Recover terminal adjoint
        adjointSketch_->reconstruct(*lhist_[0],static_cast<int>(Nt_)-2);
        // Recover terminal state sensitivity
        whist_[1]->set(*whist_[0]);
        stateSensSketch_->reconstruct(*whist_[0],static_cast<int>(Nt_)-3);
      }
      // Update dynamic constraint and objective
      con_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(Nt_-1), timeStamp_[Nt_-1]);
      obj_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(Nt_-1), timeStamp_[Nt_-1]);
      // Solve for terminal condition
      setTerminalConditionHess(*phist_[lindex],
                               *vp.get(Nt_-1),    *lhist_[lindex],
                               *whist_[uindex-1], *whist_[uindex],
                               *uhist_[uindex-1], *uhist_[uindex],
                               *xp.get(Nt_-1),    timeStamp_[Nt_-1]);
      if (useSymHess_) {
        // Add mixed derivative Hessian
        addMixedHessLag(*hvp.get(Nt_-1),   *lhist_[lindex],
                        *whist_[uindex-1], *whist_[uindex],
                        *uhist_[uindex-1], *uhist_[uindex],
                        *xp.get(Nt_-1),    timeStamp_[Nt_-1]);
      }
      // Add adjoint sensitivity to Hessian
      addAdjointSens(*hvp.get(Nt_-1),   *phist_[lindex],
                     *uhist_[uindex-1], *uhist_[uindex],
                     *xp.get(Nt_-1),    timeStamp_[Nt_-1]);
      // Run reverse time stepper
      for (size_type k = Nt_-2; k > 0; --k) {
        // Compute new components of rhs
        computeNewStateHessLag(*rhs_,             *lhist_[lindex],
                               *whist_[uindex-1], *whist_[uindex],
                               *uhist_[uindex-1], *uhist_[uindex],
                               *xp.get(k+1),      timeStamp_[k+1],
                               false);
        computeNewMixedHessLag(*rhs_,
                               *vp.get(k+1),      *lhist_[lindex],
                               *uhist_[uindex-1], *uhist_[uindex],
                               *xp.get(k+1),      timeStamp_[k+1],
                               true);
        computeNewStateJacobian(*rhs_,             *phist_[lindex],
                                *uhist_[uindex-1], *uhist_[uindex],
                                *xp.get(k+1),      timeStamp_[k+1],
                                true);
        // Recover state, adjoint and state sensitivity from sketch
        if (useSketch_) {
          uhist_[1]->set(*uhist_[0]);
          whist_[1]->set(*whist_[0]);
          if (k==1) {
            uhist_[0]->set(*u0_);
            whist_[0]->zero();
          }
          else {
            stateSketch_->reconstruct(*uhist_[0],static_cast<int>(k)-2);
            stateSensSketch_->reconstruct(*whist_[0],static_cast<int>(k)-2);
          }
          adjointSketch_->reconstruct(*lhist_[0],static_cast<int>(k)-1);
        }
        uindex = (useSketch_ ? 1 : k);
        lindex = (useSketch_ ? 0 : k);
        // Update dynamic constraint and objective
        con_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(k), timeStamp_[k]);
        obj_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(k), timeStamp_[k]);
        // Compute old components of rhs
        computeOldStateHessLag(*rhs_,             *lhist_[lindex],
                               *whist_[uindex-1], *whist_[uindex],
                               *uhist_[uindex-1], *uhist_[uindex],
                               *xp.get(k),        timeStamp_[k],
                               true);
        computeOldMixedHessLag(*rhs_,
                               *vp.get(k),        *lhist_[lindex],
                               *uhist_[uindex-1], *uhist_[uindex],
                               *xp.get(k),        timeStamp_[k],
                               true);
        // Solve for adjoint on interval k
        advanceAdjointSens(*phist_[lindex],   *rhs_,
                           *uhist_[uindex-1], *uhist_[uindex],
                           *xp.get(k),        timeStamp_[k]);
        if (useSymHess_) {
          // Add mixed derivative Hessian
          addMixedHessLag(*hvp.get(k),       *lhist_[lindex],
                          *whist_[uindex-1], *whist_[uindex],
                          *uhist_[uindex-1], *uhist_[uindex],
                          *xp.get(k),        timeStamp_[k]);
        }
        // Add adjoint sensitivity to Hessian
        addAdjointSens(*hvp.get(k),       *phist_[lindex],
                       *uhist_[uindex-1], *uhist_[uindex],
                       *xp.get(k),        timeStamp_[k]);
      }
    }
    else {
      Objective<Real>::hessVec(hv,v,x,tol);
    }
  }

private:
  /***************************************************************************/
  /************ Method to solve state equation *******************************/
  /***************************************************************************/
  void solveState(const Vector<Real> &x) {
    if (!isStateComputed_) {
      const Real one(1);
      const PartitionedVector<Real> &xp = partition(x);
      // Set initial condition
      uhist_[0]->set(*u0_);
      if (useSketch_) {
        //stateSketch_->advance(one,*uhist_[0],0,one);
      }
      // Run time stepper
      size_type index;
      for (size_type k = 1; k < Nt_; ++k) {
        index = (useSketch_ ? 1 : k);
        // Update dynamic constraint and objective
        con_->update(*uhist_[index-1], *uhist_[index], *xp.get(k), timeStamp_[k]);
        obj_->update(*uhist_[index-1], *uhist_[index], *xp.get(k), timeStamp_[k]);
        // Solve state on current time interval
        con_->solve(*cprimal_, *uhist_[index-1], *uhist_[index], *xp.get(k), timeStamp_[k]);
        // Sketch state
        if (useSketch_) {
          stateSketch_->advance(one,*uhist_[index],static_cast<int>(k)-1,one);
          uhist_[index-1]->set(*uhist_[index]);
        }
      }
      isStateComputed_ = true;
    }
  }

  Real updateSketch(const Vector<Real> &x, const Real tol) {
    const PartitionedVector<Real> &xp = partition(x);
    Real err(0), cnorm(0);
    bool flag = true;
    while (flag) {
      err = static_cast<Real>(0);
      uhist_[0]->set(*u0_);
      for (size_type k = 1; k < Nt_; ++k) {
        stateSketch_->reconstruct(*uhist_[1],static_cast<int>(k));
        con_->value(*cprimal_, *uhist_[0], *uhist_[1], *xp.get(k), timeStamp_[k]);
        cnorm = cprimal_->norm();
        err   = (cnorm > err ? cnorm : err); // Use Linf norm.  Could change to use, e.g., L2.
        if (err > tol) {
          break;
        }
      }
      if (err > tol) {
        rankState_ *= updateFactor_; // Perhaps there is a better update strategy
        rankState_  = (maxRank_ < rankState_ ? maxRank_ : rankState_);
        stateSketch_->setRank(rankState_);
        solveState(x);
      }
      else {
        flag = false;
        break;
      }
    }
    return err;
  }

  /***************************************************************************/
  /************ Methods to solve adjoint equation ****************************/
  /***************************************************************************/
  void setTerminalCondition(Vector<Real> &l,
                      const Vector<Real> &uold, const Vector<Real> &unew,
                      const Vector<Real> &z,    const TimeStamp<Real> &ts) {
    obj_->gradient_un(*udual_, uold, unew, z, ts);
    con_->applyInverseAdjointJacobian_un(l, *udual_, uold, unew, z, ts);
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

  void solveAdjoint(const Vector<Real> &x) {
    if (!isAdjointComputed_) {
      const PartitionedVector<Real> &xp = partition(x);
      const Real one(1);
      size_type uindex = (useSketch_ ? 1 : Nt_-1);
      size_type lindex = (useSketch_ ? 0 : Nt_-1);
      // Must first compute solve the state equation
      solveState(x);
      // Recover state from sketch
      if (useSketch_) {
        uhist_[1]->set(*uhist_[0]);
        stateSketch_->reconstruct(*uhist_[0],static_cast<int>(Nt_)-3);
      }
      // Update dynamic constraint and objective
      con_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(Nt_-1), timeStamp_[Nt_-1]);
      obj_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(Nt_-1), timeStamp_[Nt_-1]);
      // Solve for terminal condition
      setTerminalCondition(*lhist_[lindex],
                           *uhist_[uindex-1], *uhist_[uindex],
                           *xp.get(Nt_-1),    timeStamp_[Nt_-1]);
      if (useSketch_) {
        adjointSketch_->advance(one,*lhist_[lindex],static_cast<int>(Nt_)-2,one);
      }
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
            stateSketch_->reconstruct(*uhist_[0],static_cast<int>(k)-2);
          }
        }
        uindex = (useSketch_ ? 1 : k);
        lindex = (useSketch_ ? 0 : k);
        // Update dynamic constraint and objective
        con_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(k), timeStamp_[k]);
        obj_->update(*uhist_[uindex-1], *uhist_[uindex], *xp.get(k), timeStamp_[k]);
        // Solve for adjoint on interval k
        advanceAdjoint(*lhist_[lindex],   *rhs_,
                       *uhist_[uindex-1], *uhist_[uindex],
                       *xp.get(k),        timeStamp_[k]);
        if (useSketch_) {
          adjointSketch_->advance(one,*lhist_[lindex],static_cast<int>(k)-1,one);
        }
      }
      isAdjointComputed_ = true;
    }
  }

  /***************************************************************************/
  /************ Method for gradient computation ******************************/
  /***************************************************************************/
  void updateGradient(Vector<Real> &g,    const Vector<Real> &l,
                const Vector<Real> &uold, const Vector<Real> &unew,
                const Vector<Real> &z,    const TimeStamp<Real> &ts) {
    const Real one(1);
    obj_->gradient_z(g, uold, unew, z, ts);
    con_->applyAdjointJacobian_z(*zdual_, l, uold, unew, z, ts);
    g.axpy(-one,*zdual_);
  }

  /***************************************************************************/
  /************ Method to solve state sensitivity equation *******************/
  /***************************************************************************/
  void advanceStateSens(Vector<Real>    &wnew, const Vector<Real> &v,
                  const Vector<Real>    &wold, const Vector<Real> &uold,
                  const Vector<Real>    &unew, const Vector<Real> &z,
                  const TimeStamp<Real> &ts) {
    const Real one(1);
    con_->applyJacobian_z(*crhs_, v, uold, unew, z, ts);
    con_->applyJacobian_uo(*cprimal_, wold, uold, unew, z, ts);
    crhs_->axpy(-one, *cprimal_);
    con_->applyInverseJacobian_un(wnew, *crhs_, uold, unew, z, ts);
  }

  /***************************************************************************/
  /************ Methods to solve adjoint sensitivity equation ****************/
  /***************************************************************************/
  void setTerminalConditionHess(Vector<Real> &p,
                          const Vector<Real> &v,    const Vector<Real>    &l,
                          const Vector<Real> &wold, const Vector<Real>    &wnew,
                          const Vector<Real> &uold, const Vector<Real>    &unew,
                          const Vector<Real> &z,    const TimeStamp<Real> &ts) {
    const Real one(1);
    // Mixed derivative rhs term
    con_->applyAdjointHessian_z_un(*rhs_, l, v, uold, unew, z, ts);
    obj_->hessVec_un_z(*udual_, v, uold, unew, z, ts);
    rhs_->axpy(-one, *udual_);
    // State derivative rhs term
    con_->applyAdjointHessian_un_un(*udual_, l, wnew, uold, unew, z, ts);
    rhs_->axpy(-one, *udual_);
    obj_->hessVec_un_un(*udual_, wnew, uold, unew, z, ts);
    rhs_->plus(*udual_);
    con_->applyAdjointHessian_uo_un(*udual_, l, wold, uold, unew, z, ts);
    rhs_->axpy(-one, *udual_);
    obj_->hessVec_un_uo(*udual_, wold, uold, unew, z, ts);
    rhs_->plus(*udual_);
    // Invert adjoint Jacobian
    con_->applyInverseAdjointJacobian_un(p, *rhs_, uold, unew, z, ts);
  }

  void computeOldStateHessLag(Vector<Real> &Hv,   const Vector<Real>    &l,
                        const Vector<Real> &wold, const Vector<Real>    &wnew,
                        const Vector<Real> &uold, const Vector<Real>    &unew,
                        const Vector<Real> &z,    const TimeStamp<Real> &ts,
                        const bool sumInto = false) {
    const Real one(1);
    // Compute new/old Hessian of Lagrangian
    if (!sumInto) {
      obj_->hessVec_un_uo(Hv, wold, uold, unew, z, ts);
    }
    else {
      obj_->hessVec_un_uo(*udual_, wold, uold, unew, z, ts);
      Hv.plus(*udual_);
    }
    con_->applyAdjointHessian_uo_un(*udual_, l, wold, uold, unew, z, ts);
    Hv.axpy(-one,*udual_);
    // Compute new/new Hessian of Lagrangian
    obj_->hessVec_un_un(*udual_, wnew, uold, unew, z, ts);
    Hv.plus(*udual_);
    con_->applyAdjointHessian_un_un(*udual_, l, wnew, uold, unew, z, ts);
    Hv.axpy(-one,*udual_);
  }

  void computeOldMixedHessLag(Vector<Real> &Hv,
                        const Vector<Real> &v,    const Vector<Real>    &l,
                        const Vector<Real> &uold, const Vector<Real>    &unew,
                        const Vector<Real> &z,    const TimeStamp<Real> &ts,
                        const bool sumInto = false) {
    const Real one(1);
    // Compute new/old Hessian of Lagrangian
    if (!sumInto) {
      con_->applyAdjointHessian_z_un(Hv, l, v, uold, unew, z, ts);
    }
    else {
      con_->applyAdjointHessian_z_un(*udual_, l, v, uold, unew, z, ts);
      Hv.plus(*udual_);
    }
    obj_->hessVec_un_z(*udual_, v, uold, unew, z, ts);
    Hv.axpy(-one, *udual_);
  }

  void computeNewStateJacobian(Vector<Real> &Hv,   const Vector<Real>    &p,
                         const Vector<Real> &uold, const Vector<Real>    &unew,
                         const Vector<Real> &z,    const TimeStamp<Real> &ts,
                         const bool sumInto = false) {
    const Real one(1);
    if (!sumInto) {
      con_->applyAdjointJacobian_uo(Hv, p, uold, unew, z, ts);
      Hv.scale(-one);
    }
    else {
      con_->applyAdjointJacobian_uo(*udual_, p, uold, unew, z, ts);
      Hv.axpy(-one, *udual_);
    }
  }

  void computeNewStateHessLag(Vector<Real> &Hv,   const Vector<Real>    &l,
                        const Vector<Real> &wold, const Vector<Real>    &wnew,
                        const Vector<Real> &uold, const Vector<Real>    &unew,
                        const Vector<Real> &z,    const TimeStamp<Real> &ts,
                        const bool sumInto = false) {
    const Real one(1);
    // Compute old/new Hessian of Lagrangian
    if (!sumInto) {
      obj_->hessVec_uo_un(Hv, wnew, uold, unew, z, ts);
    }
    else {
      obj_->hessVec_uo_un(*udual_, wnew, uold, unew, z, ts);
      Hv.plus(*udual_);
    }
    con_->applyAdjointHessian_un_uo(*udual_, l, wnew, uold, unew, z, ts);
    Hv.axpy(-one,*udual_);
    // Compute old/old Hessian of Lagrangian
    obj_->hessVec_uo_uo(*udual_, wold, uold, unew, z, ts);
    Hv.plus(*udual_);
    con_->applyAdjointHessian_uo_uo(*udual_, l, wold, uold, unew, z, ts);
    Hv.axpy(-one,*udual_);
  }

  void computeNewMixedHessLag(Vector<Real> &Hv,
                        const Vector<Real> &v,    const Vector<Real>    &l,
                        const Vector<Real> &uold, const Vector<Real>    &unew,
                        const Vector<Real> &z,    const TimeStamp<Real> &ts,
                        const bool sumInto = false) {
    const Real one(1);
    // Compute new/old Hessian of Lagrangian
    if (!sumInto) {
      con_->applyAdjointHessian_z_uo(Hv, l, v, uold, unew, z, ts);
    }
    else {
      con_->applyAdjointHessian_z_uo(*udual_, l, v, uold, unew, z, ts);
      Hv.plus(*udual_);
    }
    obj_->hessVec_uo_z(*udual_, v, uold, unew, z, ts);
    Hv.axpy(-one,*udual_);
  }

  void advanceAdjointSens(Vector<Real> &p,          Vector<Real>    &rhs,
                    const Vector<Real> &uold, const Vector<Real>    &unew,
                    const Vector<Real> &z,    const TimeStamp<Real> &ts) {
    // Solve adjoint sensitivity on current time interval
    con_->applyInverseAdjointJacobian_un(p, rhs, uold, unew, z, ts);
  }

  /***************************************************************************/
  /************ Method for Hessian-times-a-vector computation ****************/
  /***************************************************************************/
  void computeControlHessLag(Vector<Real> &Hv,
                       const Vector<Real> &v,    const Vector<Real>    &l,
                       const Vector<Real> &uold, const Vector<Real>    &unew,
                       const Vector<Real> &z,    const TimeStamp<Real> &ts) {
    const Real one(1);
    // Compute Hessian of Lagrangian
    obj_->hessVec_z_z(Hv, v, uold, unew, z, ts);
    con_->applyAdjointHessian_z_z(*zdual_, l, v, uold, unew, z, ts);
    Hv.axpy(-one, *zdual_);
  }

  void addMixedHessLag(Vector<Real> &Hv,   const Vector<Real>    &l,
                 const Vector<Real> &wold, const Vector<Real>    &wnew,
                 const Vector<Real> &uold, const Vector<Real>    &unew,
                 const Vector<Real> &z,    const TimeStamp<Real> &ts) {
    const Real one(1);
    // Compute Hessian of Lagrangian on previous time interval
    obj_->hessVec_z_uo(*zdual_, wold, uold, unew, z, ts);
    Hv.axpy(-one, *zdual_);
    con_->applyAdjointHessian_uo_z(*zdual_, l, wold, uold, unew, z, ts);
    Hv.plus(*zdual_);
    // Compute Hessian of Lagrangian on current time interval
    obj_->hessVec_z_un(*zdual_, wnew, uold, unew, z, ts);
    Hv.axpy(-one, *zdual_);
    con_->applyAdjointHessian_un_z(*zdual_, l, wnew, uold, unew, z, ts);
    Hv.plus(*zdual_);
  }

  void addAdjointSens(Vector<Real> &Hv,   const Vector<Real>    &p,
                const Vector<Real> &uold, const Vector<Real>    &unew,
                const Vector<Real> &z,    const TimeStamp<Real> &ts) {
    con_->applyAdjointJacobian_z(*zdual_, p, uold, unew, z, ts);
    Hv.plus(*zdual_);
  }
};

} // namespace ROL

#endif // ROL_REDUCEDDYNAMICOBJECTIVE_HPP
