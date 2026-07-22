// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_RISK_BOUND_CONSTRAINT_H
#define ROL_RISK_BOUND_CONSTRAINT_H

#include "ROL_StdBoundConstraint.hpp"
#include "ROL_RiskVector.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template <class Real>
class RiskBoundConstraint : public BoundConstraint<Real> {
private:
  Ptr<BoundConstraint<Real>> bc_;

  Ptr<StdBoundConstraint<Real>> statObj_bc_;
  std::vector<Real> lowerObj_, upperObj_;

  std::vector<Ptr<StdBoundConstraint<Real>>> statCon_bc_;
  std::vector<std::vector<Real>> lowerCon_, upperCon_;

  bool augmentedObj_, activatedObj_;
  int nStatObj_;

  bool augmentedCon_;
  std::vector<bool> activatedCon_;
  std::vector<int> nStatCon_;

  mutable bool isLOinitialized_, isHIinitialized_;
  mutable Ptr<RiskVector<Real>> lo_, hi_;

  void setBoundInfo(ParameterList &parlist,
                    int &nStat,
                    std::vector<Real> &lower,
                    std::vector<Real> &upper,
                    bool &augmented,
                    bool &activated) {
    lower.clear(); upper.clear();
    // Get stochastic optimization information
    std::string optType = parlist.sublist("SOL").get("Type","Risk Averse");
    if ( optType == "Risk Averse" ||
         optType == "Deviation"   ||
         optType == "Regret"      ||
         optType == "Error"       ||
         optType == "Probability" ) {
      std::string name;
      RandVarFunctionalInfo<Real>(parlist,name,nStat,lower,upper,activated);
      augmented = (nStat > 0) ? true : false;
    }
    else if ( optType == "Risk Neutral" || optType == "Mean Value" ) {
      augmented = false;
      activated = false;
      nStat     = 0;
    }
    else {
      ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
        ">>> (ROL::RiskBoundConstraint): Invalid stochastic optimization type!" << optType);
    }
  }

  bool buildObjStatBnd(Ptr<ParameterList> &parlist) {
    // Objective statistic bound
    if (parlist != nullPtr) {
      setBoundInfo(*parlist,nStatObj_,lowerObj_,upperObj_,augmentedObj_,activatedObj_);
      // Build statistic bound constraint
      if ( augmentedObj_ ) {
        statObj_bc_ = makePtr<StdBoundConstraint<Real>>(lowerObj_,upperObj_);
      }
    }
    else {
      augmentedObj_ = false;
      activatedObj_ = false;
      nStatObj_     = 0;
      statObj_bc_   = nullPtr;
    }
    // Determine whether or not bound constraint is activated
    if ( !activatedObj_ ) {
      if ( statObj_bc_ != nullPtr ) {
        statObj_bc_->deactivate();
      }
    }
    return activatedObj_;
  }

  bool buildConStatBnd(std::vector<Ptr<ParameterList>> &parlist) {
    // Constraint statistic bound
    int size = parlist.size();
    nStatCon_.clear(); nStatCon_.resize(size,0);
    lowerCon_.clear(); lowerCon_.resize(size);
    upperCon_.clear(); upperCon_.resize(size);
    activatedCon_.clear(); activatedCon_.resize(size,false);
    statCon_bc_.clear(); statCon_bc_.resize(size,nullPtr);
    bool activated = false;
    for (int i = 0; i < size; ++i) {
      if ( parlist[i] != nullPtr ) {
        bool augmented = false;
        int nStat = 0;
        std::vector<Real> lo, up;
        bool act = false;
        setBoundInfo(*parlist[i],nStat,lo,up,augmented,act);
        nStatCon_[i]     = nStat;
        lowerCon_[i]     = lo;
        upperCon_[i]     = up;
        activatedCon_[i] = act;
        augmentedCon_ = (augmented ? true : augmentedCon_);
        // Build statistic bound constraint
        if ( augmented ) {
          statCon_bc_[i] = makePtr<StdBoundConstraint<Real>>(lowerCon_[i],upperCon_[i]);
        }
      }
      else {
        activatedCon_[i] = false;
        nStatCon_[i]     = 0;
        statCon_bc_[i]   = nullPtr;
      }
      if ( !activatedCon_[i] ) {
        if ( statCon_bc_[i] != nullPtr ) {
          statCon_bc_[i]->deactivate();
        }
      }
      activated = (activatedCon_[i] ? true : activated);
    }
    return activated;
  }

public:

  // Objective risk only
  RiskBoundConstraint(Ptr<ParameterList>          &parlist,
                const Ptr<BoundConstraint<Real>>  &bc = nullPtr)
   : BoundConstraint<Real>(), bc_(bc), statObj_bc_(nullPtr),
     augmentedObj_(false), activatedObj_(false),
     augmentedCon_(false),
     isLOinitialized_(false), isHIinitialized_(false) {
    bool activatedObj = buildObjStatBnd(parlist);
    // Determine whether or not bound constraint is activated
    BoundConstraint<Real>::activate();
    if ( !activatedObj ) {
      if ( bc == nullPtr || (bc != nullPtr && !bc->isActivated()) ) {
        BoundConstraint<Real>::deactivate();
      }
    }
  }

  // Constraint risk only
  RiskBoundConstraint(std::vector<Ptr<ParameterList>> &parlist,
                const Ptr<BoundConstraint<Real>>      &bc = nullPtr)
   : BoundConstraint<Real>(), bc_(bc), statObj_bc_(nullPtr),
     augmentedObj_(false), activatedObj_(false),
     augmentedCon_(false),
     isLOinitialized_(false), isHIinitialized_(false) {
    bool activatedCon = buildConStatBnd(parlist);
    // Determine whether or not bound constraint is activated
    BoundConstraint<Real>::activate();
    if ( !activatedCon ) {
      if ( bc == nullPtr || (bc != nullPtr && !bc->isActivated()) ) {
        BoundConstraint<Real>::deactivate();
      }
    }
  }

  // Objective and constraint risk
  RiskBoundConstraint(Ptr<ParameterList>              &parlistObj,
                      std::vector<Ptr<ParameterList>> &parlistCon,
                const Ptr<BoundConstraint<Real>>      &bc = nullPtr)
   : BoundConstraint<Real>(), bc_(bc), statObj_bc_(nullPtr),
     augmentedObj_(false), activatedObj_(false),
     augmentedCon_(false),
     isLOinitialized_(false), isHIinitialized_(false) {
    bool activatedObj = buildObjStatBnd(parlistObj);
    bool activatedCon = buildConStatBnd(parlistCon);
    // Determine whether or not bound constraint is activated
    BoundConstraint<Real>::activate();
    if ( !activatedObj && !activatedCon ) {
      if ( bc == nullPtr || (bc != nullPtr && !bc->isActivated()) ) {
        BoundConstraint<Real>::deactivate();
      }
    }
  }

  // Objective only -- no statistic
  RiskBoundConstraint(const Ptr<BoundConstraint<Real>> &bc)
   : BoundConstraint<Real>(), bc_(bc), statObj_bc_(nullPtr),
     augmentedObj_(false), activatedObj_(false),
     augmentedCon_(false),
     isLOinitialized_(false), isHIinitialized_(false) {
    activatedObj_ = bc_->isActivated();
    BoundConstraint<Real>::activate();
    if (!activatedObj_) {
      BoundConstraint<Real>::deactivate();
    }
  }

  void project( Vector<Real> &x ) {
    if ( augmentedObj_ && activatedObj_ ) {
      Ptr<StdVector<Real>> xs = dynamic_cast<RiskVector<Real>&>(x).getStatisticVector(0);
      statObj_bc_->project(*xs);
    }
    if (augmentedCon_) {
      int size = statCon_bc_.size();
      for (int i = 0; i < size; ++i) {
        if (activatedCon_[i]) {
          Ptr<StdVector<Real>> xs = dynamic_cast<RiskVector<Real>&>(x).getStatisticVector(1,i);
          statCon_bc_[i]->project(*xs);
        }
      }
    }
    if ( bc_ != nullPtr && bc_->isActivated() ) {
      Ptr<Vector<Real>> xvec = dynamic_cast<RiskVector<Real>&>(x).getVector();
      bc_->project(*xvec);
    }
  }

  void projectInterior( Vector<Real> &x ) {
    if ( augmentedObj_ && activatedObj_ ) {
      Ptr<StdVector<Real>> xs = dynamic_cast<RiskVector<Real>&>(x).getStatisticVector(0);
      statObj_bc_->projectInterior(*xs);
    }
    if (augmentedCon_) {
      int size = statCon_bc_.size();
      for (int i = 0; i < size; ++i) {
        if (activatedCon_[i]) {
          Ptr<StdVector<Real>> xs = dynamic_cast<RiskVector<Real>&>(x).getStatisticVector(1,i);
          statCon_bc_[i]->projectInterior(*xs);
        }
      }
    }
    if ( bc_ != nullPtr && bc_->isActivated() ) {
      Ptr<Vector<Real>> xvec = dynamic_cast<RiskVector<Real>&>(x).getVector();
      bc_->projectInterior(*xvec);
    }
  }

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0) ) {
    if ( augmentedObj_ && activatedObj_ ) {
      Ptr<StdVector<Real>>       vs = dynamic_cast<RiskVector<Real>&>(v).getStatisticVector(0);
      Ptr<const StdVector<Real>> xs = dynamic_cast<const RiskVector<Real>&>(x).getStatisticVector(0);
      statObj_bc_->pruneUpperActive(*vs,*xs,eps);
    }
    if (augmentedCon_) {
      int size = statCon_bc_.size();
      for (int i = 0; i < size; ++i) {
        if (activatedCon_[i]) {
          Ptr<StdVector<Real>>       vs = dynamic_cast<RiskVector<Real>&>(v).getStatisticVector(1,i);
          Ptr<const StdVector<Real>> xs = dynamic_cast<const RiskVector<Real>&>(x).getStatisticVector(1,i);
          statCon_bc_[i]->pruneUpperActive(*vs,*xs,eps);
        }
      }
    }
    if ( bc_ != nullPtr && bc_->isActivated() ) {
      Ptr<Vector<Real>>       vv = dynamic_cast<RiskVector<Real>&>(v).getVector();
      Ptr<const Vector<Real>> xv = dynamic_cast<const RiskVector<Real>&>(x).getVector();
      bc_->pruneUpperActive(*vv,*xv,eps);
    }
  }

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0) ) {
    if ( augmentedObj_ && activatedObj_ ) {
      Ptr<StdVector<Real>>       vs = dynamic_cast<RiskVector<Real>&>(v).getStatisticVector(0);
      Ptr<const StdVector<Real>> gs = dynamic_cast<const RiskVector<Real>&>(g).getStatisticVector(0);
      Ptr<const StdVector<Real>> xs = dynamic_cast<const RiskVector<Real>&>(x).getStatisticVector(0);
      statObj_bc_->pruneUpperActive(*vs,*gs,*xs,xeps,geps);
    }
    if (augmentedCon_) {
      int size = statCon_bc_.size();
      for (int i = 0; i < size; ++i) {
        if (activatedCon_[i]) {
          Ptr<StdVector<Real>>       vs = dynamic_cast<RiskVector<Real>&>(v).getStatisticVector(1,i);
          Ptr<const StdVector<Real>> gs = dynamic_cast<const RiskVector<Real>&>(g).getStatisticVector(1,i);
          Ptr<const StdVector<Real>> xs = dynamic_cast<const RiskVector<Real>&>(x).getStatisticVector(1,i);
          statCon_bc_[i]->pruneUpperActive(*vs,*gs,*xs,xeps,geps);
        }
      }
    }
    if ( bc_ != nullPtr && bc_->isActivated() ) {
      Ptr<Vector<Real>>       vv = dynamic_cast<RiskVector<Real>&>(v).getVector();
      Ptr<const Vector<Real>> gv = dynamic_cast<const RiskVector<Real>&>(g).getVector();
      Ptr<const Vector<Real>> xv = dynamic_cast<const RiskVector<Real>&>(x).getVector();
      bc_->pruneUpperActive(*vv,*gv,*xv,xeps,geps);
    }
  }

  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &x, Real eps = Real(0) ) {
    if ( augmentedObj_ && activatedObj_ ) {
      Ptr<StdVector<Real>>       vs = dynamic_cast<RiskVector<Real>&>(v).getStatisticVector(0);
      Ptr<const StdVector<Real>> xs = dynamic_cast<const RiskVector<Real>&>(x).getStatisticVector(0);
      statObj_bc_->pruneLowerActive(*vs,*xs,eps);
    }
    if (augmentedCon_) {
      int size = statCon_bc_.size();
      for (int i = 0; i < size; ++i) {
        if (activatedCon_[i]) {
          Ptr<StdVector<Real>>       vs = dynamic_cast<RiskVector<Real>&>(v).getStatisticVector(1,i);
          Ptr<const StdVector<Real>> xs = dynamic_cast<const RiskVector<Real>&>(x).getStatisticVector(1,i);
          statCon_bc_[i]->pruneLowerActive(*vs,*xs,eps);
        }
      }
    }
    if ( bc_ != nullPtr && bc_->isActivated() ) {
      Ptr<Vector<Real>>       vv = dynamic_cast<RiskVector<Real>&>(v).getVector();
      Ptr<const Vector<Real>> xv = dynamic_cast<const RiskVector<Real>&>(x).getVector();
      bc_->pruneLowerActive(*vv,*xv,eps);
    }
  }

  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps = Real(0), Real geps = Real(0) ) {
    if ( augmentedObj_ && activatedObj_ ) {
      Ptr<StdVector<Real>>       vs = dynamic_cast<RiskVector<Real>&>(v).getStatisticVector(0);
      Ptr<const StdVector<Real>> gs = dynamic_cast<const RiskVector<Real>&>(g).getStatisticVector(0);
      Ptr<const StdVector<Real>> xs = dynamic_cast<const RiskVector<Real>&>(x).getStatisticVector(0);
      statObj_bc_->pruneLowerActive(*vs,*gs,*xs,xeps,geps);
    }
    if (augmentedCon_) {
      int size = statCon_bc_.size();
      for (int i = 0; i < size; ++i) {
        if (activatedCon_[i]) {
          Ptr<StdVector<Real>>       vs = dynamic_cast<RiskVector<Real>&>(v).getStatisticVector(1,i);
          Ptr<const StdVector<Real>> gs = dynamic_cast<const RiskVector<Real>&>(g).getStatisticVector(1,i);
          Ptr<const StdVector<Real>> xs = dynamic_cast<const RiskVector<Real>&>(x).getStatisticVector(1,i);
          statCon_bc_[i]->pruneLowerActive(*vs,*gs,*xs,xeps,geps);
        }
      }
    }
    if ( bc_ != nullPtr && bc_->isActivated() ) {
      Ptr<Vector<Real>>       vv = dynamic_cast<RiskVector<Real>&>(v).getVector();
      Ptr<const Vector<Real>> gv = dynamic_cast<const RiskVector<Real>&>(g).getVector();
      Ptr<const Vector<Real>> xv = dynamic_cast<const RiskVector<Real>&>(x).getVector();
      bc_->pruneLowerActive(*vv,*gv,*xv,xeps,geps);
    }
  }

  const Ptr<const Vector<Real>> getLowerBound(void) const {
    if (!isLOinitialized_) {
      const Ptr<const Vector<Real>> vlo = bc_->getLowerBound();
      Ptr<std::vector<Real>> lowerObj = makePtr<std::vector<Real>>(lowerObj_);
      int size = statCon_bc_.size();
      std::vector<Ptr<std::vector<Real>>> lowerCon(size);
      for (int i = 0; i < size; ++i) {
        lowerCon[i] = makePtr<std::vector<Real>>(lowerCon_[i]);
      }
      lo_ = makePtr<RiskVector<Real>>(constPtrCast<Vector<Real>>(vlo),
                                              lowerObj,lowerCon);
      isLOinitialized_ = true;
    }
    return lo_;
  }

  const Ptr<const Vector<Real>> getUpperBound(void) const {
    if (!isHIinitialized_) {
      const Ptr<const Vector<Real>> vhi = bc_->getUpperBound();
      Ptr<std::vector<Real>> upperObj = makePtr<std::vector<Real>>(upperObj_);
      int size = statCon_bc_.size();
      std::vector<Ptr<std::vector<Real>>> upperCon(size);
      for (int i = 0; i < size; ++i) {
        upperCon[i] = makePtr<std::vector<Real>>(upperCon_[i]);
      }
      hi_ = makePtr<RiskVector<Real>>(constPtrCast<Vector<Real>>(vhi),
                                              upperObj,upperCon);
      isHIinitialized_ = true;
    }
    return hi_;
  }

  bool isFeasible( const Vector<Real> &v ) {
    bool flagstat = true, flagcon = true, flagvec = true;
    if ( augmentedObj_ && activatedObj_ ) {
      Ptr<const StdVector<Real>> vs = dynamic_cast<const RiskVector<Real>&>(v).getStatisticVector(0);
      flagstat = statObj_bc_->isFeasible(*vs);
    }
    if (augmentedCon_) {
      int size = statCon_bc_.size();
      for (int i = 0; i < size; ++i) {
        if (activatedCon_[i]) {
          Ptr<const StdVector<Real>> vs = dynamic_cast<const RiskVector<Real>&>(v).getStatisticVector(1,i);
          flagcon = (!statCon_bc_[i]->isFeasible(*vs) ? false : flagcon);
        }
      }
    }
    if ( bc_ != nullPtr && bc_->isActivated() ) {
      Ptr<const Vector<Real>> vv = dynamic_cast<const RiskVector<Real>&>(v).getVector();
      flagvec = bc_->isFeasible(*vv);
    }
    return (flagstat && flagcon && flagvec);
  }

  void applyInverseScalingFunction(Vector<Real> &dv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &g) const {
    if ( augmentedObj_ && activatedObj_ ) {
      Ptr<StdVector<Real>>      dvs = dynamic_cast<RiskVector<Real>&>(dv).getStatisticVector(0);
      Ptr<const StdVector<Real>> vs = dynamic_cast<const RiskVector<Real>&>(v).getStatisticVector(0);
      Ptr<const StdVector<Real>> xs = dynamic_cast<const RiskVector<Real>&>(x).getStatisticVector(0);
      Ptr<const StdVector<Real>> gs = dynamic_cast<const RiskVector<Real>&>(g).getStatisticVector(0);
      statObj_bc_->applyInverseScalingFunction(*dvs,*vs,*xs,*gs);
    }
    if (augmentedCon_) {
      int size = statCon_bc_.size();
      for (int i = 0; i < size; ++i) {
        if (activatedCon_[i]) {
          Ptr<StdVector<Real>>      dvs = dynamic_cast<RiskVector<Real>&>(dv).getStatisticVector(1,i);
          Ptr<const StdVector<Real>> vs = dynamic_cast<const RiskVector<Real>&>(v).getStatisticVector(1,i);
          Ptr<const StdVector<Real>> xs = dynamic_cast<const RiskVector<Real>&>(x).getStatisticVector(1,i);
          Ptr<const StdVector<Real>> gs = dynamic_cast<const RiskVector<Real>&>(g).getStatisticVector(1,i);
          statCon_bc_[i]->applyInverseScalingFunction(*dvs,*vs,*xs,*gs);
        }
      }
    }
    if ( bc_ != nullPtr && bc_->isActivated() ) {
      Ptr<Vector<Real>>      dvs = dynamic_cast<RiskVector<Real>&>(dv).getVector();
      Ptr<const Vector<Real>> vs = dynamic_cast<const RiskVector<Real>&>(v).getVector();
      Ptr<const Vector<Real>> xs = dynamic_cast<const RiskVector<Real>&>(x).getVector();
      Ptr<const Vector<Real>> gs = dynamic_cast<const RiskVector<Real>&>(g).getVector();
      bc_->applyInverseScalingFunction(*dvs,*vs,*xs,*gs);
    }
  }

  void applyScalingFunctionJacobian(Vector<Real> &dv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &g) const {
    if ( augmentedObj_ && activatedObj_ ) {
      Ptr<StdVector<Real>>      dvs = dynamic_cast<RiskVector<Real>&>(dv).getStatisticVector(0);
      Ptr<const StdVector<Real>> vs = dynamic_cast<const RiskVector<Real>&>(v).getStatisticVector(0);
      Ptr<const StdVector<Real>> xs = dynamic_cast<const RiskVector<Real>&>(x).getStatisticVector(0);
      Ptr<const StdVector<Real>> gs = dynamic_cast<const RiskVector<Real>&>(g).getStatisticVector(0);
      statObj_bc_->applyScalingFunctionJacobian(*dvs,*vs,*xs,*gs);
    }
    if (augmentedCon_) {
      int size = statCon_bc_.size();
      for (int i = 0; i < size; ++i) {
        if (activatedCon_[i]) {
          Ptr<StdVector<Real>>      dvs = dynamic_cast<RiskVector<Real>&>(dv).getStatisticVector(1,i);
          Ptr<const StdVector<Real>> vs = dynamic_cast<const RiskVector<Real>&>(v).getStatisticVector(1,i);
          Ptr<const StdVector<Real>> xs = dynamic_cast<const RiskVector<Real>&>(x).getStatisticVector(1,i);
          Ptr<const StdVector<Real>> gs = dynamic_cast<const RiskVector<Real>&>(g).getStatisticVector(1,i);
          statCon_bc_[i]->applyScalingFunctionJacobian(*dvs,*vs,*xs,*gs);
        }
      }
    }
    if ( bc_ != nullPtr && bc_->isActivated() ) {
      Ptr<Vector<Real>>      dvs = dynamic_cast<RiskVector<Real>&>(dv).getVector();
      Ptr<const Vector<Real>> vs = dynamic_cast<const RiskVector<Real>&>(v).getVector();
      Ptr<const Vector<Real>> xs = dynamic_cast<const RiskVector<Real>&>(x).getVector();
      Ptr<const Vector<Real>> gs = dynamic_cast<const RiskVector<Real>&>(g).getVector();
      bc_->applyScalingFunctionJacobian(*dvs,*vs,*xs,*gs);
    }
  }

}; // class RiskBoundConstraint

} // namespace ROL

#endif
