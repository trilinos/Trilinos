// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file
    \brief  Contains definitions for std::vector bound constraints.
    \author Created by D. Ridzal and D. Kouri.
 */

#ifndef ROL_STDBOUNDCONSTRAINT_DEF_HPP
#define ROL_STDBOUNDCONSTRAINT_DEF_HPP

namespace ROL {

template<class Real>
StdBoundConstraint<Real>::StdBoundConstraint(std::vector<Real> &x, bool isLower, Real scale, Real feasTol)
  : scale_(scale), feasTol_(feasTol), min_diff_(ROL_INF<Real>()){
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

  lower_ = makePtr<StdVector<Real>>(makePtrFromRef(x_lo_));
  upper_ = makePtr<StdVector<Real>>(makePtrFromRef(x_up_));
}

template<class Real>
StdBoundConstraint<Real>::StdBoundConstraint(std::vector<Real> &l, std::vector<Real> &u, Real scale, Real feasTol)
  : x_lo_(l), x_up_(u), scale_(scale), feasTol_(feasTol) {
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
    const Real eps = feasTol_;
    const Real tol = 100.0*ROL_EPSILON<Real>();
    const Real one(1);
    if ( BoundConstraint<Real>::isLowerActivated() ) {
      for ( int i = 0; i < dim_; ++i ) {
        Real val = ((x_lo_[i] < -tol) ? (one-eps)*x_lo_[i]
                 : ((x_lo_[i] >  tol) ? (one+eps)*x_lo_[i]
                 : x_lo_[i]+eps));
        val = std::min(x_lo_[i]+eps*min_diff_, val);
        (*ex)[i] = ((*ex)[i] < val) ? val : (*ex)[i];
      }
    }
    if ( BoundConstraint<Real>::isUpperActivated() ) {
      for ( int i = 0; i < dim_; ++i ) {
        Real val = ((x_up_[i] < -tol) ? (one+eps)*x_up_[i]
                 : ((x_up_[i] >  tol) ? (one-eps)*x_up_[i]
                 : x_up_[i]-eps));
        val = std::max(x_up_[i]-eps*min_diff_, val);
        (*ex)[i] = ((*ex)[i] > val) ? val : (*ex)[i];
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
void StdBoundConstraint<Real>::buildScalingFunction(Vector<Real> &d, const Vector<Real> &x, const Vector<Real> &g) const {
  Ptr<std::vector<Real>> ed =
    dynamic_cast<StdVector<Real>&>(d).getVector();
  Ptr<const std::vector<Real>> ex =
    dynamic_cast<const StdVector<Real>&>(x).getVector();
  Ptr<const std::vector<Real>> eg =
    dynamic_cast<const StdVector<Real>&>(g).getVector();

  Real grad, lodiff, updiff, c;

  for ( int i = 0; i < dim_; ++i ) {
    grad = (*eg)[i];
    lodiff = (*ex)[i] - x_lo_[i];
    updiff = x_up_[i] - (*ex)[i];
    c = buildC(i);
    if (-grad > lodiff) {
      if (lodiff <= updiff) {
        (*ed)[i] = std::min(std::abs(grad), c);
        continue;
      }
    }
    if (+grad > updiff) {
      if (updiff <= lodiff) {
        (*ed)[i] = std::min(std::abs(grad), c);
        continue;
      }
    }
    (*ed)[i] = std::min({lodiff, updiff, c});
  }
}

template<class Real>
void StdBoundConstraint<Real>::applyInverseScalingFunction(Vector<Real> &dv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &g) const {
  buildScalingFunction(dv, x, g);

  Ptr<std::vector<Real>> edv =
    dynamic_cast<StdVector<Real>&>(dv).getVector();
  Ptr<const std::vector<Real>> ev =
    dynamic_cast<const StdVector<Real>&>(v).getVector();

  for ( int i = 0; i < dim_; ++i ) {
    (*edv)[i] = (*ev)[i]/(*edv)[i];
  }
}

template<class Real>
void StdBoundConstraint<Real>::applyScalingFunctionJacobian(Vector<Real> &dv, const Vector<Real> &v, const Vector<Real> &x, const Vector<Real> &g) const {
  buildScalingFunction(dv, x, g);

  Ptr<std::vector<Real>> edv =
    dynamic_cast<StdVector<Real>&>(dv).getVector();
  Ptr<const std::vector<Real>> ev =
    dynamic_cast<const StdVector<Real>&>(v).getVector();
  Ptr<const std::vector<Real>> ex =
    dynamic_cast<const StdVector<Real>&>(x).getVector();
  Ptr<const std::vector<Real>> eg =
    dynamic_cast<const StdVector<Real>&>(g).getVector();

  Real zero(0), one(1), indicator, d1prime;

  for ( int i = 0; i < dim_; ++i ) {
    indicator = (*edv)[i] < buildC(i) ? one : zero;

    if (indicator == zero) {
      (*edv)[i] = zero;
      continue;
    }

    // When indicator is not zero...

    d1prime = sgn((*eg)[i]);
    if (d1prime == zero) {
      d1prime = one;
      if (x_up_[i] - (*ex)[i] < (*ex)[i] - x_lo_[i])
        d1prime = -one;
    }
    (*edv)[i] = d1prime*(*eg)[i]*(*ev)[i];
  }
}

}// End ROL Namespace

#endif
