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

#include "ROL_BoundConstraint.hpp"
#include "ROL_RiskVector.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template <class Real>
class RiskBoundConstraint : public BoundConstraint<Real> {
private:
  Teuchos::RCP<BoundConstraint<Real> > bc_;
  bool augmented_;
  bool activated_;
  int nStat_;
  std::vector<Real> lower_;
  std::vector<Real> upper_;

public:

  RiskBoundConstraint(Teuchos::ParameterList &parlist,
                const Teuchos::RCP<BoundConstraint<Real> > &bc = Teuchos::null)
   : BoundConstraint<Real>(), bc_(bc), augmented_(false), activated_(false), nStat_(0) {
    lower_.clear(); upper_.clear();
    Real zero(0);
    std::string optType = parlist.sublist("SOL").get("Stochastic Optimization Type","Risk Averse");
    if ( optType == "BPOE" ) {
      augmented_ = true;
      activated_ = true;
      nStat_     = 1;
      lower_.resize(nStat_,ROL_NINF<Real>());
      upper_.resize(nStat_,ROL_INF<Real>());
      lower_[0] = zero;
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

    if ( !activated_ ) {
      if ( bc == Teuchos::null || (bc != Teuchos::null && !bc->isActivated()) ) {
        BoundConstraint<Real>::deactivate();
      }
    }
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    if ( bc_ != Teuchos::null && bc_->isActivated() ) {
      Teuchos::RCP<const Vector<Real> > xv
        = (Teuchos::dyn_cast<RiskVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      bc_->update(*xv,flag,iter);
    }
  }

  void project( Vector<Real> &x ) {
    if ( augmented_ && activated_ ) {
      std::vector<Real> xstat(nStat_,0);
      for ( int i = 0; i < nStat_; i++ ) {
        xstat[i] = Teuchos::dyn_cast<RiskVector<Real> >(x).getStatistic(i);
        xstat[i] = std::min(upper_[i],std::max(lower_[i],xstat[i]));
      }
      (Teuchos::dyn_cast<RiskVector<Real> >(x)).setStatistic(xstat);
    }
    if ( bc_ != Teuchos::null && bc_->isActivated() ) {
      Teuchos::RCP<Vector<Real> > xvec = Teuchos::dyn_cast<RiskVector<Real> >(x).getVector();
      bc_->project(*xvec);
      (Teuchos::dyn_cast<RiskVector<Real> >(x)).setVector(*xvec);
    }
  }

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) {
    if ( augmented_ && activated_ ) {
      Real xstat(0);
      std::vector<Real> vstat(nStat_,0);
      for (int i = 0; i < nStat_; i++) {
        xstat = Teuchos::dyn_cast<const RiskVector<Real> >(x).getStatistic(i);
        if ( xstat >= upper_[i] - eps ) {
          vstat[i] = (Real)0;
        }
        else {
          vstat[i] = Teuchos::dyn_cast<RiskVector<Real> >(v).getStatistic(i);
        }
      }
      Teuchos::dyn_cast<RiskVector<Real> >(v).setStatistic(vstat);
    }
    if ( bc_ != Teuchos::null && bc_->isActivated() ) {
      Teuchos::RCP<Vector<Real> > vvec = Teuchos::dyn_cast<RiskVector<Real> >(v).getVector();
      Teuchos::RCP<const Vector<Real> > xvec = Teuchos::dyn_cast<const RiskVector<Real> >(x).getVector();
      bc_->pruneUpperActive(*vvec,*xvec,eps);
    }
  }

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) {
    if ( augmented_ && activated_ ) {
      Real gstat(0), xstat(0);
      std::vector<Real> vstat(nStat_,0);
      for (int i = 0; i < nStat_; i++) {
        gstat = Teuchos::dyn_cast<const RiskVector<Real> >(g).getStatistic(i);
        xstat = Teuchos::dyn_cast<const RiskVector<Real> >(x).getStatistic(i);
        if ( (xstat >= upper_[i] - eps) && (gstat < (Real)0) ) {
          vstat[i] = (Real)0;
        }
        else {
          vstat[i] = Teuchos::dyn_cast<RiskVector<Real> >(v).getStatistic(i);
        }
      }
      Teuchos::dyn_cast<RiskVector<Real> >(v).setStatistic(vstat);
    }
    if ( bc_ != Teuchos::null && bc_->isActivated() ) {
      Teuchos::RCP<Vector<Real> > vvec = Teuchos::dyn_cast<RiskVector<Real> >(v).getVector();
      Teuchos::RCP<const Vector<Real> > gvec = Teuchos::dyn_cast<const RiskVector<Real> >(g).getVector();
      Teuchos::RCP<const Vector<Real> > xvec = Teuchos::dyn_cast<const RiskVector<Real> >(x).getVector();
      bc_->pruneUpperActive(*vvec,*gvec,*xvec,eps);
    }
  }
 
  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) {
    if ( augmented_ && activated_ ) {
      Real xstat(0);
      std::vector<Real> vstat(nStat_,0);
      for (int i = 0; i < nStat_; i++) {
        xstat = Teuchos::dyn_cast<const RiskVector<Real> >(x).getStatistic(i);
        if ( xstat <= lower_[i] + eps ) {
          vstat[i] = (Real)0;
        }
        else {
          vstat[i] = Teuchos::dyn_cast<RiskVector<Real> >(v).getStatistic(i);
        }
      }
      Teuchos::dyn_cast<RiskVector<Real> >(v).setStatistic(vstat);
    }
    if ( bc_ != Teuchos::null && bc_->isActivated() ) {
      Teuchos::RCP<Vector<Real> > vvec = Teuchos::dyn_cast<RiskVector<Real> >(v).getVector();
      Teuchos::RCP<const Vector<Real> > xvec = Teuchos::dyn_cast<const RiskVector<Real> >(x).getVector();
      bc_->pruneLowerActive(*vvec,*xvec,eps);
    }
  }

  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) {
    if ( augmented_ && activated_ ) {
      Real gstat(0), xstat(0);
      std::vector<Real> vstat(nStat_,0);
      for (int i = 0; i < nStat_; i++) {
        gstat = Teuchos::dyn_cast<const RiskVector<Real> >(g).getStatistic(i);
        xstat = Teuchos::dyn_cast<const RiskVector<Real> >(x).getStatistic(i);
        if ( (xstat <= lower_[i] + eps) && (gstat > (Real)0) ) {
          vstat[i] = (Real)0;
        }
        else {
          vstat[i] = Teuchos::dyn_cast<RiskVector<Real> >(v).getStatistic(i);
        }
      }
      Teuchos::dyn_cast<RiskVector<Real> >(v).setStatistic(vstat);
    }
    if ( bc_ != Teuchos::null && bc_->isActivated() ) {
      Teuchos::RCP<Vector<Real> > vvec = Teuchos::dyn_cast<RiskVector<Real> >(v).getVector();
      Teuchos::RCP<const Vector<Real> > gvec = Teuchos::dyn_cast<const RiskVector<Real> >(g).getVector();
      Teuchos::RCP<const Vector<Real> > xvec = Teuchos::dyn_cast<const RiskVector<Real> >(x).getVector();
      bc_->pruneLowerActive(*vvec,*gvec,*xvec,eps);
    }
  } 

  void setVectorToUpperBound( Vector<Real> &u ) {
    if ( augmented_ && activated_ ) {
      Teuchos::dyn_cast<RiskVector<Real> >(u).setStatistic(upper_);
    }
    if ( bc_ != Teuchos::null && bc_->isActivated() ) {
      Teuchos::RCP<Vector<Real> > uvec = Teuchos::dyn_cast<RiskVector<Real> >(u).getVector();
      bc_->setVectorToUpperBound(*uvec);
    }
  }

  void setVectorToLowerBound( Vector<Real> &l ) {
    if ( augmented_ && activated_ ) {
      Teuchos::dyn_cast<RiskVector<Real> >(l).setStatistic(lower_);
    }
    if ( bc_ != Teuchos::null && bc_->isActivated() ) {
      Teuchos::RCP<Vector<Real> > lvec = Teuchos::dyn_cast<RiskVector<Real> >(l).getVector();
      bc_->setVectorToLowerBound(*lvec);
    }
  }

  void pruneActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) {
    if ( augmented_ && activated_ ) {
      Real xstat(0);
      std::vector<Real> vstat(nStat_,0);
      for (int i = 0; i < nStat_; i++) {
        xstat = Teuchos::dyn_cast<const RiskVector<Real> >(x).getStatistic(i);
        if ( (xstat <= lower_[i] + eps) || (xstat >= upper_[i] - eps) ) {
          vstat[i] = (Real)0;
        }
        else {
          vstat[i] = Teuchos::dyn_cast<RiskVector<Real> >(v).getStatistic(i);
        }
      }
      Teuchos::dyn_cast<RiskVector<Real> >(v).setStatistic(vstat);
    }
    if ( bc_ != Teuchos::null && bc_->isActivated() ) {
      Teuchos::RCP<Vector<Real> > vvec = Teuchos::dyn_cast<RiskVector<Real> >(v).getVector();
      Teuchos::RCP<const Vector<Real> > xvec = Teuchos::dyn_cast<const RiskVector<Real> >(x).getVector();
      bc_->pruneActive(*vvec,*xvec,eps);
    }
  }

  void pruneActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) {
    if ( augmented_ && activated_ ) {
      Real gstat(0), xstat(0);
      std::vector<Real> vstat(nStat_,0);
      for (int i = 0; i < nStat_; i++) {
        gstat = Teuchos::dyn_cast<const RiskVector<Real> >(g).getStatistic(i);
        xstat = Teuchos::dyn_cast<const RiskVector<Real> >(x).getStatistic(i);
        if ( ((xstat <= lower_[i] + eps) && (gstat > (Real)0)) ||
             ((xstat >= upper_[i] - eps) && (gstat < (Real)0)) ) {
          vstat[i] = (Real)0;
        }
        else {
          vstat[i] = Teuchos::dyn_cast<RiskVector<Real> >(v).getStatistic(i);
        }
      }
      Teuchos::dyn_cast<RiskVector<Real> >(v).setStatistic(vstat);
    }
    if ( bc_ != Teuchos::null && bc_->isActivated() ) {
      Teuchos::RCP<Vector<Real> > vvec = Teuchos::dyn_cast<RiskVector<Real> >(v).getVector();
      Teuchos::RCP<const Vector<Real> > gvec = Teuchos::dyn_cast<const RiskVector<Real> >(g).getVector();
      Teuchos::RCP<const Vector<Real> > xvec = Teuchos::dyn_cast<const RiskVector<Real> >(x).getVector();
      bc_->pruneActive(*vvec,*gvec,*xvec,eps);
    }
  }

  bool isFeasible( const Vector<Real> &v ) { 
    bool flagstat = true, flagvec = true;
    if ( augmented_ && activated_ ) {
      Real vstat(0);
      for ( int i = 0; i < nStat_; i++ ) {
        vstat = Teuchos::dyn_cast<const RiskVector<Real> >(v).getStatistic(i);
        flagstat *= ((vstat >= lower_[i] && vstat <= upper_[i]) ? true : false);
      }
    }
    if ( bc_ != Teuchos::null && bc_->isActivated() ) {
      Teuchos::RCP<const Vector<Real> > vvec
        = (Teuchos::dyn_cast<RiskVector<Real> >(const_cast<Vector<Real> &>(v))).getVector();
      if ( bc_->isActivated() ) {
        flagvec = bc_->isFeasible(*vvec);
      }
    }
    return (flagstat && flagvec);
  }

}; // class RiskBoundConstraint

} // namespace ROL

#endif
