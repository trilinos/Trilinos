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
  Teuchos::RCP<BoundConstraint<Real> > bc_;
  Teuchos::RCP<StdBoundConstraint<Real> > stat_bc_;
  std::vector<Real> lower_, upper_;

  bool augmented_, activated_;
  int nStat_;

  mutable bool isLOinitialized_, isHIinitialized_; 
  mutable Teuchos::RCP<RiskVector<Real> > lo_, hi_;

public:

  RiskBoundConstraint(Teuchos::ParameterList &parlist,
                const Teuchos::RCP<BoundConstraint<Real> > &bc = Teuchos::null)
   : BoundConstraint<Real>(), bc_(bc), stat_bc_(Teuchos::null),
     augmented_(false), activated_(false), nStat_(0),
     isLOinitialized_(false), isHIinitialized_(false) {
    lower_.clear(); upper_.clear();
    // Get stochastic optimization information
    std::string optType = parlist.sublist("SOL").get("Stochastic Optimization Type","Risk Averse");
    if ( optType == "BPOE" ) {
      augmented_ = true;
      activated_ = true;
      nStat_     = 1;
      lower_.resize(nStat_,ROL_NINF<Real>());
      upper_.resize(nStat_,ROL_INF<Real>());
      lower_[0] = static_cast<Real>(0);
    }
    else if ( optType == "Risk Averse" ) {
      std::string name;
      RiskMeasureInfo<Real>(parlist,name,nStat_,lower_,upper_,activated_);
      augmented_ = (nStat_ > 0) ? true : false;
    }
    else if ( optType == "Risk Neutral" || optType == "Mean Value" ) {
      augmented_ = false;
      activated_ = false;
      nStat_     = 0;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
        ">>> (ROL::RiskBoundConstraint): Invalid stochastic optimization type!" << optType);
    }
    // Build statistic bound constraint
    if ( augmented_ ) {
      stat_bc_ = Teuchos::rcp(new StdBoundConstraint<Real>(lower_,upper_));
    }
    // Determine whether or not bound constraint is activated
    BoundConstraint<Real>::activate();
    if ( !activated_ ) {
      if ( stat_bc_ != Teuchos::null ) {
        stat_bc_->deactivate();
      }
      if ( bc == Teuchos::null || (bc != Teuchos::null && !bc->isActivated()) ) {
        BoundConstraint<Real>::deactivate();
      }
    }
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    if ( augmented_ && activated_ ) {
      Teuchos::RCP<const StdVector<Real> > xs = Teuchos::dyn_cast<const RiskVector<Real> >(x).getStatistic();
      stat_bc_->update(*xs,flag,iter);
    }
    if ( bc_ != Teuchos::null && bc_->isActivated() ) {
      Teuchos::RCP<const Vector<Real> > xv = Teuchos::dyn_cast<const RiskVector<Real> >(x).getVector();
      bc_->update(*xv,flag,iter);
    }
  }

  void project( Vector<Real> &x ) {
    if ( augmented_ && activated_ ) {
      Teuchos::RCP<StdVector<Real> > xs = Teuchos::dyn_cast<RiskVector<Real> >(x).getStatistic();
      stat_bc_->project(*xs);
    }
    if ( bc_ != Teuchos::null && bc_->isActivated() ) {
      Teuchos::RCP<Vector<Real> > xvec = Teuchos::dyn_cast<RiskVector<Real> >(x).getVector();
      bc_->project(*xvec);
    }
  }

  void projectInterior( Vector<Real> &x ) {
    if ( augmented_ && activated_ ) {
      Teuchos::RCP<StdVector<Real> > xs = Teuchos::dyn_cast<RiskVector<Real> >(x).getStatistic();
      stat_bc_->projectInterior(*xs);
    }
    if ( bc_ != Teuchos::null && bc_->isActivated() ) {
      Teuchos::RCP<Vector<Real> > xvec = Teuchos::dyn_cast<RiskVector<Real> >(x).getVector();
      bc_->projectInterior(*xvec);
    }
  }

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0 ) {
    if ( augmented_ && activated_ ) {
      Teuchos::RCP<StdVector<Real> >       vs = Teuchos::dyn_cast<RiskVector<Real> >(v).getStatistic();
      Teuchos::RCP<const StdVector<Real> > xs = Teuchos::dyn_cast<const RiskVector<Real> >(x).getStatistic();
      stat_bc_->pruneUpperActive(*vs,*xs,eps);
    }
    if ( bc_ != Teuchos::null && bc_->isActivated() ) {
      Teuchos::RCP<Vector<Real> >       vv = Teuchos::dyn_cast<RiskVector<Real> >(v).getVector();
      Teuchos::RCP<const Vector<Real> > xv = Teuchos::dyn_cast<const RiskVector<Real> >(x).getVector();
      bc_->pruneUpperActive(*vv,*xv,eps);
    }
  }

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0 ) {
    if ( augmented_ && activated_ ) {
      Teuchos::RCP<StdVector<Real> >       vs = Teuchos::dyn_cast<RiskVector<Real> >(v).getStatistic();
      Teuchos::RCP<const StdVector<Real> > gs = Teuchos::dyn_cast<const RiskVector<Real> >(g).getStatistic();
      Teuchos::RCP<const StdVector<Real> > xs = Teuchos::dyn_cast<const RiskVector<Real> >(x).getStatistic();
      stat_bc_->pruneUpperActive(*vs,*gs,*xs,eps);
    }
    if ( bc_ != Teuchos::null && bc_->isActivated() ) {
      Teuchos::RCP<Vector<Real> >       vv = Teuchos::dyn_cast<RiskVector<Real> >(v).getVector();
      Teuchos::RCP<const Vector<Real> > gv = Teuchos::dyn_cast<const RiskVector<Real> >(g).getVector();
      Teuchos::RCP<const Vector<Real> > xv = Teuchos::dyn_cast<const RiskVector<Real> >(x).getVector();
      bc_->pruneUpperActive(*vv,*gv,*xv,eps);
    }
  }
 
  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0 ) {
    if ( augmented_ && activated_ ) {
      Teuchos::RCP<StdVector<Real> >       vs = Teuchos::dyn_cast<RiskVector<Real> >(v).getStatistic();
      Teuchos::RCP<const StdVector<Real> > xs = Teuchos::dyn_cast<const RiskVector<Real> >(x).getStatistic();
      stat_bc_->pruneLowerActive(*vs,*xs,eps);
    }
    if ( bc_ != Teuchos::null && bc_->isActivated() ) {
      Teuchos::RCP<Vector<Real> >       vv = Teuchos::dyn_cast<RiskVector<Real> >(v).getVector();
      Teuchos::RCP<const Vector<Real> > xv = Teuchos::dyn_cast<const RiskVector<Real> >(x).getVector();
      bc_->pruneLowerActive(*vv,*xv,eps);
    }
  }

  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0 ) {
    if ( augmented_ && activated_ ) {
      Teuchos::RCP<StdVector<Real> >       vs = Teuchos::dyn_cast<RiskVector<Real> >(v).getStatistic();
      Teuchos::RCP<const StdVector<Real> > gs = Teuchos::dyn_cast<const RiskVector<Real> >(g).getStatistic();
      Teuchos::RCP<const StdVector<Real> > xs = Teuchos::dyn_cast<const RiskVector<Real> >(x).getStatistic();
      stat_bc_->pruneLowerActive(*vs,*gs,*xs,eps);
    }
    if ( bc_ != Teuchos::null && bc_->isActivated() ) {
      Teuchos::RCP<Vector<Real> >       vv = Teuchos::dyn_cast<RiskVector<Real> >(v).getVector();
      Teuchos::RCP<const Vector<Real> > gv = Teuchos::dyn_cast<const RiskVector<Real> >(g).getVector();
      Teuchos::RCP<const Vector<Real> > xv = Teuchos::dyn_cast<const RiskVector<Real> >(x).getVector();
      bc_->pruneLowerActive(*vv,*gv,*xv,eps);
    }
  } 

  const Teuchos::RCP<const Vector<Real> > getLowerBound(void) const {
    if (!isLOinitialized_) {
      const Teuchos::RCP<const Vector<Real> > vlo = bc_->getLowerBound();
      lo_ = Teuchos::rcp(new RiskVector<Real>(Teuchos::rcp_const_cast<Vector<Real> >(vlo),
                                              lower_, augmented_));
      isLOinitialized_ = true;
    }
    return lo_;
  }

  const Teuchos::RCP<const Vector<Real> > getUpperBound(void) const {
    if (!isHIinitialized_) {
      const Teuchos::RCP<const Vector<Real> > vhi = bc_->getUpperBound();
      hi_ = Teuchos::rcp(new RiskVector<Real>(Teuchos::rcp_const_cast<Vector<Real> >(vhi),
                                              upper_, augmented_));
      isHIinitialized_ = true;
    }
    return hi_;
  }

  bool isFeasible( const Vector<Real> &v ) { 
    bool flagstat = true, flagvec = true;
    if ( augmented_ && activated_ ) {
      Teuchos::RCP<const StdVector<Real> > vs = Teuchos::dyn_cast<const RiskVector<Real> >(v).getStatistic();
      flagstat = stat_bc_->isFeasible(*vs);
    }
    if ( bc_ != Teuchos::null && bc_->isActivated() ) {
      Teuchos::RCP<const Vector<Real> > vv = Teuchos::dyn_cast<const RiskVector<Real> >(v).getVector();
      flagvec = bc_->isFeasible(*vv);
    }
    return (flagstat && flagvec);
  }

}; // class RiskBoundConstraint

} // namespace ROL

#endif
