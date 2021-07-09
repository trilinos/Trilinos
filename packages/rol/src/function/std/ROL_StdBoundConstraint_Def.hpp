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

/** \file
    \brief  Contains definitions for std::vector bound constraints.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef ROL_STDBOUNDCONSTRAINT_DEF_HPP
#define ROL_STDBOUNDCONSTRAINT_DEF_HPP

namespace ROL {

template<class Real>
StdBoundConstraint<Real>::StdBoundConstraint(std::vector<Real> &x, bool isLower, Real scale)
  : scale_(scale) {
  dim_ = x.size();
  x_lo_.clear(); x_up_.clear();
  if (isLower) {
    x_lo_.assign(x.begin(),x.end());
    x_up_.resize(dim_,ROL_INF<Real>());
    BoundConstraint<Real>::activateLower();
  }
  else {
    x_lo_.resize(dim_,ROL_NINF<Real>());
    x_up_.assign(x.begin(),x.end());
    BoundConstraint<Real>::activateUpper();
  }
  min_diff_ = ROL_INF<Real>();

  lower_ = makePtr<StdVector<Real>>(makePtrFromRef(x_lo_));
  upper_ = makePtr<StdVector<Real>>(makePtrFromRef(x_up_));
}

template<class Real>
StdBoundConstraint<Real>::StdBoundConstraint(std::vector<Real> &l, std::vector<Real> &u, Real scale)
  : x_lo_(l), x_up_(u), scale_(scale) {
  BoundConstraint<Real>::activate();
  dim_ = x_lo_.size();
  for ( int i = 0; i < dim_; i++ ) {
    if ( i == 0 ) {
      min_diff_ = x_up_[i] - x_lo_[i];
    }
    else {
      min_diff_ = ( (min_diff_ < (x_up_[i] - x_lo_[i])) ? min_diff_ : (x_up_[i] - x_lo_[i]) );
    }
  }
  min_diff_ *= 0.5;

  lower_ = makePtr<StdVector<Real>>(makePtrFromRef(x_lo_));
  upper_ = makePtr<StdVector<Real>>(makePtrFromRef(x_up_));
}

template<class Real>
bool StdBoundConstraint<Real>::isFeasible( const Vector<Real> &x ) {
  bool lflag = true, uflag = true;
  if ( BoundConstraint<Real>::isActivated() ) {
    Ptr<const std::vector<Real>> ex =
      dynamic_cast<const StdVector<Real>&>(x).getVector();
    if ( BoundConstraint<Real>::isLowerActivated() ) {
      for ( int i = 0; i < dim_; ++i ) {
        if ( (*ex)[i] < x_lo_[i] ) {
          lflag = false;
          break;
        }
      }
    }
    if ( BoundConstraint<Real>::isUpperActivated() ) {
      for ( int i = 0; i < dim_; ++i ) {
        if ( (*ex)[i] > x_up_[i] ) {
          uflag = false;
          break;
        }
      }
    }
  }
  return (lflag && uflag);
}

template<class Real>
void StdBoundConstraint<Real>::project( Vector<Real> &x ) {
  if ( BoundConstraint<Real>::isActivated() ) {
    Ptr<std::vector<Real>> ex =
      dynamic_cast<StdVector<Real>&>(x).getVector();
    if ( BoundConstraint<Real>::isLowerActivated() ) {
      for ( int i = 0; i < dim_; ++i ) {
        (*ex)[i] = std::max(x_lo_[i],(*ex)[i]);
      }
    }
    if ( BoundConstraint<Real>::isUpperActivated() ) {
      for ( int i = 0; i < dim_; ++i ) {
        (*ex)[i] = std::min(x_up_[i],(*ex)[i]);
      }
    }
  }
}

template<class Real>
void StdBoundConstraint<Real>::projectInterior( Vector<Real> &x ) {
  if ( BoundConstraint<Real>::isActivated() ) {
    Ptr<std::vector<Real>> ex =
        dynamic_cast<StdVector<Real>&>(x).getVector();
    const Real eps(1e-1), tol(100.0*ROL_EPSILON<Real>()), one(1);
    if ( BoundConstraint<Real>::isLowerActivated() ) {
      for ( int i = 0; i < dim_; ++i ) {
        Real val = ((x_lo_[i] < -tol) ? (one-eps)*x_lo_[i]
                 : ((x_lo_[i] >  tol) ? (one+eps)*x_lo_[i]
                 : x_lo_[i]+eps));
        val = std::min(x_lo_[i]+eps*min_diff_, val);
        (*ex)[i] = ((*ex)[i] < x_lo_[i]+tol) ? val : (*ex)[i];
      }
    }
    if ( BoundConstraint<Real>::isUpperActivated() ) {
      for ( int i = 0; i < dim_; ++i ) {
        Real val = ((x_up_[i] < -tol) ? (one+eps)*x_up_[i]
                 : ((x_up_[i] >  tol) ? (one-eps)*x_up_[i]
                 : x_up_[i]-eps));
        val = std::max(x_up_[i]-eps*min_diff_, val);
        (*ex)[i] = ((*ex)[i] > x_up_[i]-tol) ? val : (*ex)[i];
      }
    }
  }
}

template<class Real>
void StdBoundConstraint<Real>::pruneLowerActive(Vector<Real> &v, const Vector<Real> &x, Real eps) {
  if ( BoundConstraint<Real>::isLowerActivated() ) {
    Ptr<const std::vector<Real>> ex =
      dynamic_cast<const StdVector<Real>&>(x).getVector();
    Ptr<std::vector<Real>> ev =
      dynamic_cast<StdVector<Real>&>(v).getVector();
    Real epsn = std::min(scale_*eps,min_diff_);
    for ( int i = 0; i < dim_; ++i ) {
      if ( ((*ex)[i] <= x_lo_[i]+epsn) ) {
        (*ev)[i] = static_cast<Real>(0);
      }
    }
  }
}

template<class Real>
void StdBoundConstraint<Real>::pruneUpperActive(Vector<Real> &v, const Vector<Real> &x, Real eps) {
  if ( BoundConstraint<Real>::isUpperActivated() ) {
    Ptr<const std::vector<Real>> ex =
      dynamic_cast<const StdVector<Real>&>(x).getVector();
    Ptr<std::vector<Real>> ev =
      dynamic_cast<StdVector<Real>&>(v).getVector();
    Real epsn = std::min(scale_*eps,min_diff_);
    for ( int i = 0; i < dim_; ++i ) {
      if ( ((*ex)[i] >= x_up_[i]-epsn) ) {
        (*ev)[i] = static_cast<Real>(0);
      }
    }
  }
}

template<class Real>
void StdBoundConstraint<Real>::pruneLowerActive(Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps, Real geps) {
  if ( BoundConstraint<Real>::isLowerActivated() ) {
    Ptr<const std::vector<Real>> ex =
      dynamic_cast<const StdVector<Real>&>(x).getVector();
    Ptr<const std::vector<Real>> eg =
      dynamic_cast<const StdVector<Real>&>(g).getVector();
    Ptr<std::vector<Real>> ev =
      dynamic_cast<StdVector<Real>&>(v).getVector();
    Real epsn = std::min(scale_*xeps,this->min_diff_);
    for ( int i = 0; i < dim_; ++i ) {
      if ( (*ex)[i] <= x_lo_[i]+epsn && (*eg)[i] > geps ) {
        (*ev)[i] = static_cast<Real>(0);
      }
    }
  }
}

template<class Real>
void StdBoundConstraint<Real>::pruneUpperActive(Vector<Real> &v, const Vector<Real> &g, const Vector<Real> &x, Real xeps, Real geps) {
  if ( BoundConstraint<Real>::isUpperActivated() ) {
    Ptr<const std::vector<Real>> ex = 
      dynamic_cast<const StdVector<Real>&>(x).getVector();
    Ptr<const std::vector<Real>> eg =
      dynamic_cast<const StdVector<Real>&>(g).getVector();
    Ptr<std::vector<Real>> ev =
      dynamic_cast<StdVector<Real>&>(v).getVector();
    Real epsn = std::min(scale_*xeps,min_diff_);
    for ( int i = 0; i < dim_; ++i ) {
      if ( (*ex)[i] >= x_up_[i]-epsn && (*eg)[i] < -geps ) {
        (*ev)[i] = static_cast<Real>(0);
      }
    }
  }
}

}// End ROL Namespace

#endif
