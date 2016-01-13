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
  Real lower_;
  Real upper_;

public:

  RiskBoundConstraint(Teuchos::ParameterList &parlist,
                const Teuchos::RCP<BoundConstraint<Real> > &bc = Teuchos::null)
   : BoundConstraint<Real>(), bc_(bc), augmented_(false), lower_(ROL_NINF), upper_(ROL_INF) {
    std::string type = parlist.sublist("SOL").sublist("Risk Measure").get("Name","CVaR");
    if ( type == "CVaR" || type == "HMCR" ||
         type == "Log-Exponential Quadrangle" ||
         type == "Quantile-Based Quadrangle" ||
         type == "Truncated Mean Quadrangle" ) {
      augmented_ = true;
    }
    if ( !(bc_->isActivated()) ) {
      BoundConstraint<Real>::deactivate();
    }
  }

  RiskBoundConstraint(const Teuchos::RCP<BoundConstraint<Real> > &bc = Teuchos::null,
                      const bool augmented = false,
                      const Real lower = ROL_NINF,
                      const Real upper = ROL_INF)
    : bc_(bc), augmented_(augmented) {
    lower_ = std::min(lower,upper);
    upper_ = std::max(lower,upper);
    if (!augmented_ && !(bc_->isActivated()) ) {
      BoundConstraint<Real>::deactivate();
    }
  }

  RiskBoundConstraint(const std::string name,
                      const Teuchos::RCP<BoundConstraint<Real> > &bc = Teuchos::null)
    : bc_(bc), augmented_(true), lower_(ROL_NINF), upper_(ROL_INF) {
    if ( name == "BPOE" ) { lower_ = 0.; }
  }

  void update( const Vector<Real> &x, bool flag = true, int iter = -1 ) {
    if ( bc_ != Teuchos::null ) {
      Teuchos::RCP<const Vector<Real> > xv
        = (Teuchos::dyn_cast<RiskVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      bc_->update(*xv,flag,iter);
    }
  }

  void project( Vector<Real> &x ) {
    if ( augmented_ ) {
      Real xvar = Teuchos::dyn_cast<RiskVector<Real> >(x).getStatistic();
      xvar = std::min(upper_,std::max(lower_,xvar));
      (Teuchos::dyn_cast<RiskVector<Real> >(x)).setStatistic(xvar);
    }
    if ( bc_ != Teuchos::null ) {
      Teuchos::RCP<Vector<Real> > xvec = Teuchos::rcp_const_cast<Vector<Real> >(
        (Teuchos::dyn_cast<RiskVector<Real> >(x)).getVector());
      bc_->project(*xvec);
      (Teuchos::dyn_cast<RiskVector<Real> >(x)).setVector(*xvec);
    }
  }

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) {
    if ( augmented_ ) {
      Real xvar = Teuchos::dyn_cast<RiskVector<Real> >(const_cast<Vector<Real> &>(x)).getStatistic();
      if ( xvar >= upper_ - eps ) {
        (Teuchos::dyn_cast<RiskVector<Real> >(v)).setStatistic(0.0);
      }
    }
    if ( bc_ != Teuchos::null ) {
      Teuchos::RCP<Vector<Real> > vvec = Teuchos::rcp_const_cast<Vector<Real> >(
        (Teuchos::dyn_cast<RiskVector<Real> >(v)).getVector());
      Teuchos::RCP<const Vector<Real> > xvec =
          (Teuchos::dyn_cast<RiskVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      bc_->pruneUpperActive(*vvec,*xvec,eps);
    }
  }

  void pruneUpperActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) {
    if ( augmented_ ) {
      Real gvar = Teuchos::dyn_cast<RiskVector<Real> >(const_cast<Vector<Real> &>(g)).getStatistic();
      Real xvar = Teuchos::dyn_cast<RiskVector<Real> >(const_cast<Vector<Real> &>(x)).getStatistic();
      if ( (xvar >= upper_ - eps) && gvar < 0.0 ) {
        (Teuchos::dyn_cast<RiskVector<Real> >(v)).setStatistic(0.0);
      }
    }
    if ( bc_ != Teuchos::null ) {
      Teuchos::RCP<Vector<Real> > vvec = Teuchos::rcp_const_cast<Vector<Real> >(
        (Teuchos::dyn_cast<RiskVector<Real> >(v)).getVector());
      Teuchos::RCP<const Vector<Real> > gvec =
          (Teuchos::dyn_cast<RiskVector<Real> >(const_cast<Vector<Real> &>(g))).getVector();
      Teuchos::RCP<const Vector<Real> > xvec =
          (Teuchos::dyn_cast<RiskVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      bc_->pruneUpperActive(*vvec,*gvec,*xvec,eps);
    }
  }
 
  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) {
    if ( augmented_ ) {
      Real xvar = Teuchos::dyn_cast<RiskVector<Real> >(const_cast<Vector<Real> &>(x)).getStatistic();
      if ( xvar <= lower_ + eps ) {
        (Teuchos::dyn_cast<RiskVector<Real> >(v)).setStatistic(0.0);
      }
    }
    if ( bc_ != Teuchos::null ) {
      Teuchos::RCP<Vector<Real> > vvec = Teuchos::rcp_const_cast<Vector<Real> >(
        (Teuchos::dyn_cast<RiskVector<Real> >(v)).getVector());
      Teuchos::RCP<const Vector<Real> > xvec
        = (Teuchos::dyn_cast<RiskVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      bc_->pruneLowerActive(*vvec,*xvec,eps);
    }
  }

  void pruneLowerActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) {
    if ( augmented_ ) {
      Real gvar = Teuchos::dyn_cast<RiskVector<Real> >(const_cast<Vector<Real> &>(g)).getStatistic();
      Real xvar = Teuchos::dyn_cast<RiskVector<Real> >(const_cast<Vector<Real> &>(x)).getStatistic();
      if ( (xvar <= lower_ + eps) && gvar > 0.0 ) {
        (Teuchos::dyn_cast<RiskVector<Real> >(v)).setStatistic(0.0);
      }
    }
    if ( bc_ != Teuchos::null ) {
      Teuchos::RCP<Vector<Real> > vvec = Teuchos::rcp_const_cast<Vector<Real> >(
        (Teuchos::dyn_cast<RiskVector<Real> >(v)).getVector());
      Teuchos::RCP<const Vector<Real> > gvec
        = (Teuchos::dyn_cast<RiskVector<Real> >(const_cast<Vector<Real> &>(g))).getVector();
      Teuchos::RCP<const Vector<Real> > xvec
        = (Teuchos::dyn_cast<RiskVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      bc_->pruneLowerActive(*vvec,*gvec,*xvec,eps);
    }
  } 

  void setVectorToUpperBound( Vector<Real> &u ) {
    if ( augmented_ ) {
      (Teuchos::dyn_cast<RiskVector<Real> >(u)).setStatistic(upper_);
    }
    if ( bc_ != Teuchos::null ) {
      Teuchos::RCP<Vector<Real> > uvec = Teuchos::rcp_const_cast<Vector<Real> >(
        (Teuchos::dyn_cast<RiskVector<Real> >(u)).getVector());
      bc_->setVectorToUpperBound(*uvec);
    }
  }

  void setVectorToLowerBound( Vector<Real> &l ) {
    if ( augmented_ ) {
      (Teuchos::dyn_cast<RiskVector<Real> >(l)).setStatistic(lower_);
    }
    if ( bc_ != Teuchos::null ) {
      Teuchos::RCP<Vector<Real> > lvec = Teuchos::rcp_const_cast<Vector<Real> >(
        (Teuchos::dyn_cast<RiskVector<Real> >(l)).getVector());
      bc_->setVectorToLowerBound(*lvec);
    }
  }

  void pruneActive( Vector<Real> &v, const Vector<Real> &x, Real eps = 0.0 ) {
    if ( augmented_ ) {
      Real xvar = Teuchos::dyn_cast<RiskVector<Real> >(const_cast<Vector<Real> &>(x)).getStatistic();
      if ( (xvar <= lower_ + eps) || (xvar >= upper_ - eps) ) {
        (Teuchos::dyn_cast<RiskVector<Real> >(v)).setStatistic(0.0);
      }
    }
    if ( bc_ != Teuchos::null ) {
      Teuchos::RCP<Vector<Real> > vvec = Teuchos::rcp_const_cast<Vector<Real> >(
        (Teuchos::dyn_cast<RiskVector<Real> >(v)).getVector());
      Teuchos::RCP<const Vector<Real> > xvec =
          (Teuchos::dyn_cast<RiskVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      bc_->pruneActive(*vvec,*xvec,eps);
    }
  }

  void pruneActive( Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real eps = 0.0 ) {
    if (augmented_) {
      Real gvar = Teuchos::dyn_cast<RiskVector<Real> >(const_cast<Vector<Real> &>(g)).getStatistic();
      Real xvar = Teuchos::dyn_cast<RiskVector<Real> >(const_cast<Vector<Real> &>(x)).getStatistic();
      if ( ((xvar <= lower_ + eps) && gvar > 0.0) ||
           ((xvar >= upper_ - eps) && gvar < 0.0) ) {
        (Teuchos::dyn_cast<RiskVector<Real> >(v)).setStatistic(0.0);
      }
    }
    if ( bc_ != Teuchos::null ) {
      Teuchos::RCP<Vector<Real> > vvec = Teuchos::rcp_const_cast<Vector<Real> >(
        (Teuchos::dyn_cast<RiskVector<Real> >(v)).getVector());
      Teuchos::RCP<const Vector<Real> > gvec =
          (Teuchos::dyn_cast<RiskVector<Real> >(const_cast<Vector<Real> &>(g))).getVector();
      Teuchos::RCP<const Vector<Real> > xvec =
          (Teuchos::dyn_cast<RiskVector<Real> >(const_cast<Vector<Real> &>(x))).getVector();
      bc_->pruneActive(*vvec,*gvec,*xvec,eps);
    }
  }

  bool isFeasible( const Vector<Real> &v ) { 
    bool flagstat = true, flagvec = true;
    if ( augmented_ ) {
      Real vvar = Teuchos::dyn_cast<RiskVector<Real> >(const_cast<Vector<Real> &>(v)).getStatistic();
      flagstat = ((vvar >= lower_ && vvar <= upper_) ? true : false);
    }
    if ( bc_ != Teuchos::null ) {
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
