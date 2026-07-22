// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SPARSEGRIDGENERATOR_DEF_HPP
#define ROL_SPARSEGRIDGENERATOR_DEF_HPP

namespace ROL {

template <class Real>
SparseGridGenerator<Real>::SparseGridGenerator(const ROL::Ptr<BatchManager<Real> > &bman, 
                                               const QuadratureInfo &info, const bool adaptive)
  : SampleGenerator<Real>(bman), info_(info), isVectorInit_(false) {
  adaptive_ = info.adaptive;
  if ( adaptive_ ) {
    indices_   = ROL::makePtr<SparseGridIndexSet<Real>>(info.dim, info.maxLevel);
    grid_      = ROL::makePtr<Quadrature<Real>>(info.dim);
    error_     = static_cast<Real>(0);
    direction_ = 0;
  }
  else {
    grid_ = ROL::makePtr<Quadrature<Real>>(info);
  }
  // Split points and weights across processors
  setSamples(true);
}

template <class Real>
SparseGridGenerator<Real>::SparseGridGenerator(const ROL::Ptr<BatchManager<Real> > &bman, 
                                               const char* SGinfo, const char* SGdata, 
                                               const bool isNormalized)
    : SampleGenerator<Real>(bman), adaptive_(false), isVectorInit_(false) {
  grid_ = ROL::makePtr<Quadrature<Real>>(SGinfo,SGdata,isNormalized);
  // Split points and weights across processors
  setSamples(true);
}

template<class Real> 
void SparseGridGenerator<Real>::buildDiffRule(Quadrature<Real> &outRule,
                                              const std::vector<int> &index) const {
  Real one(1);
  int numPoints = 0;
  for ( int i = 0; i < info_.dim; i++ ) {
    numPoints = growthRule1D(index[i],(info_.growth1D)[i],(info_.rule1D)[i]);
    Quadrature<Real> diffRule((info_.rule1D)[i],numPoints,(info_.normalized));
    if ( numPoints != 1 ) {
      numPoints = growthRule1D(index[i]-1,(info_.growth1D)[i],(info_.rule1D)[i]);
      Quadrature<Real> rule((info_.rule1D)[i],numPoints,(info_.normalized));
      diffRule.update(-one,rule);
    }
    outRule = kron_prod<Real>(outRule,diffRule);
  }
}

template<class Real>
void SparseGridGenerator<Real>::update(const Vector<Real> &x) {
  SampleGenerator<Real>::update(x);
  if ( adaptive_ ) {
    indices_->reset();
    grid_ = ROL::makePtr<Quadrature<Real>>(info_.dim);
    index_.clear();
    direction_ = 0;
    error_ = static_cast<Real>(0);
    // Split points and weights across processors
    setSamples(true);
  }
}

template<class Real>
void SparseGridGenerator<Real>::refine(void) {
  if ( adaptive_ ) {
    npts_ = 0;
    ROL::Ptr<Quadrature<Real> > rule;
//    bool terminate = false;
//int cnt = 0;
//    while (!terminate) {
      // Select index to investigate
      if ( indices_->isEmpty() ) {
        // Start from the vector of ones
        index_.resize(info_.dim,1);
        search_index_ = index_;
        direction_    = info_.dim;
      }
      else {
        if (direction_ < (info_.dim-1)) {
          // Continue investigating forward neighbors
          direction_++;
        }
        else {
          // Select index corresponding to largest error
          Real error(0);
          indices_->get(error,index_);
          error_    -= error;
          direction_ = 0;
        }
        search_index_ = index_;
        search_index_[direction_]++;
      }
      rule = ROL::makePtr<Quadrature<Real>>(info_.dim);
      if (    !(indices_->isMember(search_index_))           // Check if index is old/active
           && !(indices_->isMaxLevelExceeded(search_index_)) // Check if index violates maxLevel
           && indices_->isAdmissible(search_index_) ) {      // Check if index is admissible
        // Build difference rule
        buildDiffRule(*rule,search_index_);
        npts_ = rule->getNumPoints();
//        terminate = true;
      }
      else if ( indices_->isActiveEmpty() ) {
        npts_ = 0;
//        terminate = true;
      }
////if (info_.print) {
////std::cout << "IN REFINE: CNT = " << cnt << "  ERROR = " << error_ << std::endl;
////}
////cnt++;
//    }
    // Set values of difference rule as points and weights
    updateSamples(*rule);
  }
}

template<class Real>
Real SparseGridGenerator<Real>::computeError(std::vector<Real> &vals){
  if ( adaptive_ ) {
    if (npts_>0) {
      Real myerror(0), error(0);
      for ( int i = 0; i < SampleGenerator<Real>::numMySamples(); ++i ) {
        myerror += vals[i]*SampleGenerator<Real>::getMyWeight(i);
      }
      SampleGenerator<Real>::sumAll(&myerror,&error,1);
      error = std::abs(error);
      // Update global error and index sets
      error_ += error;
      indices_->add(error, search_index_);
    }
    // Return
    return error_;
  }
  else {
    return static_cast<Real>(0);
  }
}

template<class Real>
Real SparseGridGenerator<Real>::computeError(std::vector<ROL::Ptr<Vector<Real> > > &vals,
                                             const Vector<Real> &x ){
  if ( adaptive_ ) {
    if ( !isVectorInit_ ) {
      mydiff_ = x.dual().clone();
      diff_   = x.dual().clone();
      isVectorInit_ = true;
    }
    if (npts_>0) {
      mydiff_->zero();
      for ( int i = 0; i < SampleGenerator<Real>::numMySamples(); ++i ) {
        mydiff_->axpy(SampleGenerator<Real>::getMyWeight(i),(*vals[i]));
      }
      SampleGenerator<Real>::sumAll(*mydiff_,*diff_);
      Real error = diff_->norm();
      // Update global error and index sets
      error_ += error;
      indices_->add(error, search_index_);
    }
    // Return
    return error_;
  }
  else {
    return static_cast<Real>(0);
  }
}

template<class Real> 
void SparseGridGenerator<Real>::splitSamples(std::vector<std::vector<Real> > &mypts,
                                             std::vector<Real> &mywts) {
  // Get global points and weights
  std::vector<std::vector<Real> > pts;
  std::vector<Real> wts;
  grid_->getCubature(pts,wts);
  // Split global points and weights across processors
  int frac = pts.size()/SampleGenerator<Real>::numBatches();
  int rem  = pts.size()%SampleGenerator<Real>::numBatches();
  int npts = frac;
  if ( SampleGenerator<Real>::batchID() < rem ) {
    npts++;
  } 
  mypts.resize(npts);
  mywts.resize(npts);
  int index = 0;
  for ( int i = 0; i < npts; i++ ) {
    index = i*(SampleGenerator<Real>::numBatches()) + (SampleGenerator<Real>::batchID());
    mywts[i] = wts[index];
    mypts[i] = pts[index];
  }
}

template<class Real>
void SparseGridGenerator<Real>::setSamples(bool inConstructor) {
  if (info_.print && !inConstructor) {
    std::cout << info_.name << ": Number of Sparse-Grid Points: "
              << grid_->getNumPoints() << std::endl;
  }
  if ( adaptive_ || inConstructor ) {
    // Split samples based on PID
    std::vector<std::vector<Real> > mypts;
    std::vector<Real> mywts;
    splitSamples(mypts,mywts);
    // Set local points and weights
    SampleGenerator<Real>::setPoints(mypts);
    SampleGenerator<Real>::setWeights(mywts);
  }
}

template<class Real>
void SparseGridGenerator<Real>::updateSamples( Quadrature<Real> &grid ) {
  Real one(1);
  // Add increment to stored sparse grid
  grid_->update(one,grid);
  // Split global stored points and weights across processors
  std::vector<std::vector<Real> > mypts;
  std::vector<Real> mywts;
  splitSamples(mypts,mywts);
  // Get global incremental points and weights
  std::vector<std::vector<Real> > pts_inc;
  std::vector<Real> wts_inc;
  grid.getCubature(pts_inc,wts_inc);
  // STL set intesection of global incremental points and local stored points
  std::vector<std::vector<Real> > pts_inter;
  typename std::vector<std::vector<Real> >::iterator mypts_it;
  typename std::vector<std::vector<Real> >::iterator pts_inc_it;
  for ( pts_inc_it=pts_inc.begin(); pts_inc_it!=pts_inc.end(); pts_inc_it++ ) {
    for ( mypts_it=mypts.begin(); mypts_it!=mypts.end(); mypts_it++ ) {
      if ( *mypts_it == *pts_inc_it ) {
        pts_inter.push_back(*mypts_it);
        break;
      }
    }
  }
  // Get intersection weights
  std::vector<Real> wts_inter;
  for ( unsigned i = 0; i < pts_inter.size(); i++ ) {
    wts_inter.push_back(grid.getWeight(pts_inter[i]));
  }
  // Set local points and weights
  SampleGenerator<Real>::setPoints(pts_inter);
  SampleGenerator<Real>::setWeights(wts_inter);
}

template<class Real>
void SparseGridGenerator<Real>::printIndexSet(void) const {
  if (indices_ != ROL::nullPtr) {
    indices_->print(info_.name,SampleGenerator<Real>::batchID());
  }
}

} // namespace ROL

#endif
