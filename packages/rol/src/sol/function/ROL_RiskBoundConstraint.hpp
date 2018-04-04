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

#ifndef ROL_RISK_BOUND_CONSTRAINT_H
#define ROL_RISK_BOUND_CONSTRAINT_H

#include "ROL_StdBoundConstraint.hpp"
#include "ROL_RiskVector.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template <class Real>
class RiskBoundConstraint : public BoundConstraint<Real> {
private:
  ROL::Ptr<BoundConstraint<Real> > bc_;

  ROL::Ptr<StdBoundConstraint<Real> > statObj_bc_;
  std::vector<Real> lowerObj_, upperObj_;

  std::vector<ROL::Ptr<StdBoundConstraint<Real> > > statCon_bc_;
  std::vector<std::vector<Real> > lowerCon_, upperCon_;

  bool augmentedObj_, activatedObj_;
  int nStatObj_;

  bool augmentedCon_;
  std::vector<bool> activatedCon_;
  std::vector<int> nStatCon_;

  mutable bool isLOinitialized_, isHIinitialized_; 
  mutable ROL::Ptr<RiskVector<Real> > lo_, hi_;

  void setBoundInfo(Teuchos::ParameterList &parlist,
                    int &nStat,
                    std::vector<Real> &lower,
                    std::vector<Real> &upper,
                    bool &augmented,
                    bool &activated) {
    lower.clear(); upper.clear();
    // Get stochastic optimization information
    std::string optType = parlist.sublist("SOL").get("Stochastic Component Type","Risk Averse");
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
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
        ">>> (ROL::RiskBoundConstraint): Invalid stochastic optimization type!" << optType);
    }
  }

  bool buildObjStatBnd(Teuchos::RCP<Teuchos::ParameterList> &parlist) {
    // Objective statistic bound
    if (parlist != ROL::nullPtr) {
      setBoundInfo(*parlist,nStatObj_,lowerObj_,upperObj_,augmentedObj_,activatedObj_);
      // Build statistic bound constraint
      if ( augmentedObj_ ) {
        statObj_bc_ = ROL::makePtr<StdBoundConstraint<Real>>(lowerObj_,upperObj_);
      }
    }
    else {
      augmentedObj_ = false;
      activatedObj_ = false;
      nStatObj_     = 0;
      statObj_bc_   = ROL::nullPtr;
    }
    // Determine whether or not bound constraint is activated
    if ( !activatedObj_ ) {
      if ( statObj_bc_ != ROL::nullPtr ) {
        statObj_bc_->deactivate();
      }
    }
    return activatedObj_;
  }

  bool buildConStatBnd(std::vector<Teuchos::RCP<Teuchos::ParameterList> > &parlist) {
    // Constraint statistic bound
    int size = parlist.size();
    nStatCon_.clear(); nStatCon_.resize(size,0);
    lowerCon_.clear(); lowerCon_.resize(size);
    upperCon_.clear(); upperCon_.resize(size);
    activatedCon_.clear(); activatedCon_.resize(size,false);
    statCon_bc_.clear(); statCon_bc_.resize(size,ROL::nullPtr);
    bool activated = false;
    for (int i = 0; i < size; ++i) {
      if ( parlist[i] != ROL::nullPtr ) {
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
          statCon_bc_[i] = ROL::makePtr<StdBoundConstraint<Real>>(lowerCon_[i],upperCon_[i]);
        }
      }
      else {
        activatedCon_[i] = false;
        nStatCon_[i]     = 0;
        statCon_bc_[i]   = ROL::nullPtr;
      }
      if ( !activatedCon_[i] ) {
        if ( statCon_bc_[i] != ROL::nullPtr ) {
          statCon_bc_[i]->deactivate();
        }
      }
      activated = (activatedCon_[i] ? true : activated);
    }
    return activated;
  }

public:

  // Objective risk only
  RiskBoundConstraint(ROL::Ptr<Teuchos::ParameterList > &parlist,
                const ROL::Ptr<BoundConstraint<Real> >  &bc = ROL::nullPtr)
   : BoundConstraint<Real>(), bc_(bc), statObj_bc_(ROL::nullPtr),
     augmentedObj_(false), activatedObj_(false),
     augmentedCon_(false),
     isLOinitialized_(false), isHIinitialized_(false) {
    bool activatedObj = buildObjStatBnd(parlist);
    // Determine whether or not bound constraint is activated
    BoundConstraint<Real>::activate();
    if ( !activatedObj ) {
      if ( bc == ROL::nullPtr || (bc != ROL::nullPtr && !bc->isActivated()) ) {
        BoundConstraint<Real>::deactivate();
      }
    }
  }

  // Constraint risk only
  RiskBoundConstraint(std::vector<Teuchos::RCP<Teuchos::ParameterList> > &parlist,
                const ROL::Ptr<BoundConstraint<Real> >               &bc = ROL::nullPtr)
   : BoundConstraint<Real>(), bc_(bc), statObj_bc_(ROL::nullPtr),
     augmentedObj_(false), activatedObj_(false),
     augmentedCon_(false),
     isLOinitialized_(false), isHIinitialized_(false) {
    bool activatedCon = buildConStatBnd(parlist);
    // Determine whether or not bound constraint is activated
    BoundConstraint<Real>::activate();
    if ( !activatedCon ) {
      if ( bc == ROL::nullPtr || (bc != ROL::nullPtr && !bc->isActivated()) ) {
        BoundConstraint<Real>::deactivate();
      }
    }
  }

  // Objective and constraint risk
  RiskBoundConstraint(Teuchos::RCP<Teuchos::ParameterList>               &parlistObj,
                      std::vector<Teuchos::RCP<Teuchos::ParameterList> > &parlistCon,
                const ROL::Ptr<BoundConstraint<Real> >               &bc = ROL::nullPtr)
   : BoundConstraint<Real>(), bc_(bc), statObj_bc_(ROL::nullPtr),
     augmentedObj_(false), activatedObj_(false),
     augmentedCon_(false),
     isLOinitialized_(false), isHIinitialized_(false) {
    bool activatedObj = buildObjStatBnd(parlistObj);
    bool activatedCon = buildConStatBnd(parlistCon);
    // Determine whether or not bound constraint is activated
    BoundConstraint<Real>::activate();
    if ( !activatedObj && !activatedCon ) {
      if ( bc == ROL::nullPtr || (bc != ROL::nullPtr && !bc->isActivated()) ) {
        BoundConstraint<Real>::deactivate();
      }
    }
  }

  // Objective only -- no statistic
  RiskBoundConstraint(const Ptr<BoundConstraint<Real> > &bc)
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

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    if ( augmentedObj_ && activatedObj_ ) {
      ROL::Ptr<const StdVector<Real> > xs = dynamic_cast<const RiskVector<Real>&>(x).getStatisticVector(0);
      statObj_bc_->update(*xs,flag,iter);
    }
    if (augmentedCon_) {
      int size = statCon_bc_.size();
      for (int i = 0; i < size; ++i) {
        if (activatedCon_[i]) {
          ROL::Ptr<const StdVector<Real> > xs = dynamic_cast<const RiskVector<Real>&>(x).getStatisticVector(1,i);
          statCon_bc_[i]->update(*xs,flag,iter);
        }
      }
    }
    if ( bc_ != ROL::nullPtr && bc_->isActivated() ) {
      ROL::Ptr<const Vector<Real> > xv = dynamic_cast<const RiskVector<Real>&>(x).getVector();
      bc_->update(*xv,flag,iter);
    }
  }

  void project( Vector<Real> &x ) {
    if ( augmentedObj_ && activatedObj_ ) {
      ROL::Ptr<StdVector<Real> > xs = dynamic_cast<RiskVector<Real>&>(x).getStatisticVector(0);
      statObj_bc_->project(*xs);
    }
    if (augmentedCon_) {
      int size = statCon_bc_.size();
      for (int i = 0; i < size; ++i) {
        if (activatedCon_[i]) {
          ROL::Ptr<StdVector<Real> > xs = dynamic_cast<RiskVector<Real>&>(x).getStatisticVector(1,i);
          statCon_bc_[i]->project(*xs);
        }
      }
    }
    if ( bc_ != ROL::nullPtr && bc_->isActivated() ) {
      ROL::Ptr<Vector<Real> > xvec = dynamic_cast<RiskVector<Real>&>(x).getVector();
      bc_->project(*xvec);
    }
  }

  void projectInterior( Vector<Real> &x ) {
    if ( augmentedObj_ && activatedObj_ ) {
      ROL::Ptr<StdVector<Real> > xs = dynamic_cast<RiskVector<Real>&>(x).getStatisticVector(0);
      statObj_bc_->projectInterior(*xs);
    }
    if (augmentedCon_) {
      int size = statCon_bc_.size();
      for (int i = 0; i < size; ++i) {
        if (activatedCon_[i]) {
          ROL::Ptr<StdVector<Real> > xs = dynamic_cast<RiskVector<Real>&>(x).getStatisticVector(1,i);
          statCon_bc_[i]->projectInterior(*xs);
        }
      }
    }
    if ( bc_ != ROL::nullPtr && bc_->isActivated() ) {
      ROL::Ptr<Vector<Real> > xvec = dynamic_cast<RiskVector<Real>&>(x).getVector();
      bc_->projectInterior(*xvec);
    }
  }

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0 ) {
    if ( augmentedObj_ && activatedObj_ ) {
      ROL::Ptr<StdVector<Real> >       vs = dynamic_cast<RiskVector<Real>&>(v).getStatisticVector(0);
      ROL::Ptr<const StdVector<Real> > xs = dynamic_cast<const RiskVector<Real>&>(x).getStatisticVector(0);
      statObj_bc_->pruneUpperActive(*vs,*xs,eps);
    }
    if (augmentedCon_) {
      int size = statCon_bc_.size();
      for (int i = 0; i < size; ++i) {
        if (activatedCon_[i]) {
          ROL::Ptr<StdVector<Real> >       vs = dynamic_cast<RiskVector<Real>&>(v).getStatisticVector(1,i);
          ROL::Ptr<const StdVector<Real> > xs = dynamic_cast<const RiskVector<Real>&>(x).getStatisticVector(1,i);
          statCon_bc_[i]->pruneUpperActive(*vs,*xs,eps);
        }
      }
    }
    if ( bc_ != ROL::nullPtr && bc_->isActivated() ) {
      ROL::Ptr<Vector<Real> >       vv = dynamic_cast<RiskVector<Real>&>(v).getVector();
      ROL::Ptr<const Vector<Real> > xv = dynamic_cast<const RiskVector<Real>&>(x).getVector();
      bc_->pruneUpperActive(*vv,*xv,eps);
    }
  }

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0 ) {
    if ( augmentedObj_ && activatedObj_ ) {
      ROL::Ptr<StdVector<Real> >       vs = dynamic_cast<RiskVector<Real>&>(v).getStatisticVector(0);
      ROL::Ptr<const StdVector<Real> > gs = dynamic_cast<const RiskVector<Real>&>(g).getStatisticVector(0);
      ROL::Ptr<const StdVector<Real> > xs = dynamic_cast<const RiskVector<Real>&>(x).getStatisticVector(0);
      statObj_bc_->pruneUpperActive(*vs,*gs,*xs,eps);
    }
    if (augmentedCon_) {
      int size = statCon_bc_.size();
      for (int i = 0; i < size; ++i) {
        if (activatedCon_[i]) {
          ROL::Ptr<StdVector<Real> >       vs = dynamic_cast<RiskVector<Real>&>(v).getStatisticVector(1,i);
          ROL::Ptr<const StdVector<Real> > gs = dynamic_cast<const RiskVector<Real>&>(g).getStatisticVector(1,i);
          ROL::Ptr<const StdVector<Real> > xs = dynamic_cast<const RiskVector<Real>&>(x).getStatisticVector(1,i);
          statCon_bc_[i]->pruneUpperActive(*vs,*gs,*xs,eps);
        }
      }
    }
    if ( bc_ != ROL::nullPtr && bc_->isActivated() ) {
      ROL::Ptr<Vector<Real> >       vv = dynamic_cast<RiskVector<Real>&>(v).getVector();
      ROL::Ptr<const Vector<Real> > gv = dynamic_cast<const RiskVector<Real>&>(g).getVector();
      ROL::Ptr<const Vector<Real> > xv = dynamic_cast<const RiskVector<Real>&>(x).getVector();
      bc_->pruneUpperActive(*vv,*gv,*xv,eps);
    }
  }
 
  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0 ) {
    if ( augmentedObj_ && activatedObj_ ) {
      ROL::Ptr<StdVector<Real> >       vs = dynamic_cast<RiskVector<Real>&>(v).getStatisticVector(0);
      ROL::Ptr<const StdVector<Real> > xs = dynamic_cast<const RiskVector<Real>&>(x).getStatisticVector(0);
      statObj_bc_->pruneLowerActive(*vs,*xs,eps);
    }
    if (augmentedCon_) {
      int size = statCon_bc_.size();
      for (int i = 0; i < size; ++i) {
        if (activatedCon_[i]) {
          ROL::Ptr<StdVector<Real> >       vs = dynamic_cast<RiskVector<Real>&>(v).getStatisticVector(1,i);
          ROL::Ptr<const StdVector<Real> > xs = dynamic_cast<const RiskVector<Real>&>(x).getStatisticVector(1,i);
          statCon_bc_[i]->pruneLowerActive(*vs,*xs,eps);
        }
      }
    }
    if ( bc_ != ROL::nullPtr && bc_->isActivated() ) {
      ROL::Ptr<Vector<Real> >       vv = dynamic_cast<RiskVector<Real>&>(v).getVector();
      ROL::Ptr<const Vector<Real> > xv = dynamic_cast<const RiskVector<Real>&>(x).getVector();
      bc_->pruneLowerActive(*vv,*xv,eps);
    }
  }

  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0 ) {
    if ( augmentedObj_ && activatedObj_ ) {
      ROL::Ptr<StdVector<Real> >       vs = dynamic_cast<RiskVector<Real>&>(v).getStatisticVector(0);
      ROL::Ptr<const StdVector<Real> > gs = dynamic_cast<const RiskVector<Real>&>(g).getStatisticVector(0);
      ROL::Ptr<const StdVector<Real> > xs = dynamic_cast<const RiskVector<Real>&>(x).getStatisticVector(0);
      statObj_bc_->pruneLowerActive(*vs,*gs,*xs,eps);
    }
    if (augmentedCon_) {
      int size = statCon_bc_.size();
      for (int i = 0; i < size; ++i) {
        if (activatedCon_[i]) {
          ROL::Ptr<StdVector<Real> >       vs = dynamic_cast<RiskVector<Real>&>(v).getStatisticVector(1,i);
          ROL::Ptr<const StdVector<Real> > gs = dynamic_cast<const RiskVector<Real>&>(g).getStatisticVector(1,i);
          ROL::Ptr<const StdVector<Real> > xs = dynamic_cast<const RiskVector<Real>&>(x).getStatisticVector(1,i);
          statCon_bc_[i]->pruneLowerActive(*vs,*gs,*xs,eps);
        }
      }
    }
    if ( bc_ != ROL::nullPtr && bc_->isActivated() ) {
      ROL::Ptr<Vector<Real> >       vv = dynamic_cast<RiskVector<Real>&>(v).getVector();
      ROL::Ptr<const Vector<Real> > gv = dynamic_cast<const RiskVector<Real>&>(g).getVector();
      ROL::Ptr<const Vector<Real> > xv = dynamic_cast<const RiskVector<Real>&>(x).getVector();
      bc_->pruneLowerActive(*vv,*gv,*xv,eps);
    }
  } 

  const ROL::Ptr<const Vector<Real> > getLowerBound(void) const {
    if (!isLOinitialized_) {
      const ROL::Ptr<const Vector<Real> > vlo = bc_->getLowerBound();
      ROL::Ptr<std::vector<Real> > lowerObj = ROL::makePtr<std::vector<Real>>(lowerObj_);
      int size = statCon_bc_.size();
      std::vector<ROL::Ptr<std::vector<Real> > > lowerCon(size);
      for (int i = 0; i < size; ++i) {
        lowerCon[i] = ROL::makePtr<std::vector<Real>>(lowerCon_[i]);
      }
      lo_ = ROL::makePtr<RiskVector<Real>>(ROL::constPtrCast<Vector<Real>>(vlo),
                                              lowerObj,lowerCon);
      isLOinitialized_ = true;
    }
    return lo_;
  }

  const ROL::Ptr<const Vector<Real> > getUpperBound(void) const {
    if (!isHIinitialized_) {
      const ROL::Ptr<const Vector<Real> > vhi = bc_->getUpperBound();
      ROL::Ptr<std::vector<Real> > upperObj = ROL::makePtr<std::vector<Real>>(upperObj_);
      int size = statCon_bc_.size();
      std::vector<ROL::Ptr<std::vector<Real> > > upperCon(size);
      for (int i = 0; i < size; ++i) {
        upperCon[i] = ROL::makePtr<std::vector<Real>>(upperCon_[i]);
      }
      hi_ = ROL::makePtr<RiskVector<Real>>(ROL::constPtrCast<Vector<Real>>(vhi),
                                              upperObj,upperCon);
      isHIinitialized_ = true;
    }
    return hi_;
  }

  bool isFeasible( const Vector<Real> &v ) { 
    bool flagstat = true, flagcon = true, flagvec = true;
    if ( augmentedObj_ && activatedObj_ ) {
      ROL::Ptr<const StdVector<Real> > vs = dynamic_cast<const RiskVector<Real>&>(v).getStatisticVector(0);
      flagstat = statObj_bc_->isFeasible(*vs);
    }
    if (augmentedCon_) {
      int size = statCon_bc_.size();
      for (int i = 0; i < size; ++i) {
        if (activatedCon_[i]) {
          ROL::Ptr<const StdVector<Real> > vs = dynamic_cast<const RiskVector<Real>&>(v).getStatisticVector(1,i);
          flagcon = (!statCon_bc_[i]->isFeasible(*vs) ? false : flagcon);
        }
      }
    }
    if ( bc_ != ROL::nullPtr && bc_->isActivated() ) {
      ROL::Ptr<const Vector<Real> > vv = dynamic_cast<const RiskVector<Real>&>(v).getVector();
      flagvec = bc_->isFeasible(*vv);
    }
    return (flagstat && flagcon && flagvec);
  }

}; // class RiskBoundConstraint

} // namespace ROL

#endif
