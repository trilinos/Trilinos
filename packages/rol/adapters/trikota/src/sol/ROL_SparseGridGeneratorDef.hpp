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

#ifndef ROL_SPARSEGRIDGENERATOR_DEF_HPP
#define ROL_SPARSEGRIDGENERATOR_DEF_HPP

namespace ROL {

template <class Real>
SparseGridGenerator<Real>::SparseGridGenerator(const Teuchos::RCP<BatchManager<Real> > &bman, 
                                               const QuadratureInfo &info, const bool adaptive)
  : SampleGenerator<Real>(bman), info_(info) {
  adaptive_ = info.adaptive;
  if ( adaptive_ ) {
    Real zero(0);
    direction_ = 0;
    error_ = zero;
    index_.clear(); activeIndex_.clear(); oldIndex_.clear();
    //index_.resize(info.dim,1);
    grid_ = Teuchos::rcp(new Quadrature<Real>(info.dim));
  }
  else {
    grid_ = Teuchos::rcp(new Quadrature<Real>(info));
    // Split points and weights across processors
    setSamples(true);
  }
}

template <class Real>
SparseGridGenerator<Real>::SparseGridGenerator(const Teuchos::RCP<BatchManager<Real> > &bman, 
                                               const char* SGinfo, const char* SGdata, 
                                               const bool isNormalized)
    : SampleGenerator<Real>(bman), adaptive_(false) {
  grid_ = Teuchos::rcp(new Quadrature<Real>(SGinfo,SGdata,isNormalized));
  // Split points and weights across processors
  setSamples(true);
}

template<class Real>
bool SparseGridGenerator<Real>::checkMaxLevel(std::vector<int> &index) {
  int  level = 0;
  bool useMax = true;
  for ( unsigned l = 0; l < index.size(); l++ ) {
    if ( useMax ) {
      level = std::max(level,index[l]);
    }
    else {
      level += index[l];
    }
  }
  if ( level >= info_.maxLevel ) {
    return true;
  }
  return false;
}

template<class Real>
bool SparseGridGenerator<Real>::isAdmissible(std::vector<int> &index, int direction) {
  for ( int i = 0; i < info_.dim; i++ ) {
    if ( index[i] > 1 && i != direction ) {
      index[i]--;
      if ( !(oldIndex_.count(index)) ) {
        return false;
      }
      index[i]++;
    }
  }
  return true;
}

template<class Real> 
void SparseGridGenerator<Real>::buildDiffRule(Quadrature<Real> &outRule,
                                              std::vector<int> &index) {
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
    Real zero(0);
    //index_.resize(info_.dim,1);
    index_.clear();
    direction_ = 0;
    error_ = zero;
    activeIndex_.clear();
    oldIndex_.clear();
    grid_ = Teuchos::rcp(new Quadrature<Real>(info_.dim));
  }
}

template<class Real>
Real SparseGridGenerator<Real>::computeError(std::vector<Real> &vals){
  if ( adaptive_ ) {
    Real myerror(0);
    for ( int i = 0; i < SampleGenerator<Real>::numMySamples(); i++ ) {
      myerror += vals[i]*SampleGenerator<Real>::getMyWeight(i);
    }
    Real error(0);
    SampleGenerator<Real>::sumAll(&myerror,&error,1);
    error = std::abs(error);
    // Update global error and index sets
    error_ += error;
    if ( activeIndex_.end() != activeIndex_.begin() ) {
      activeIndex_.insert(activeIndex_.end()--,
        std::pair<Real,std::vector<int> >(error,search_index_));
    }
    else {
      activeIndex_.insert(std::pair<Real,std::vector<int> >(error,search_index_));
    }
    // Return
    return error_;
  }
  else {
    return static_cast<Real>(0);
  }
}

template<class Real>
Real SparseGridGenerator<Real>::computeError(std::vector<Teuchos::RCP<Vector<Real> > > &vals,
                                             const Vector<Real> &x ){
  if ( adaptive_ ) {
    Teuchos::RCP<Vector<Real> > mydiff = x.clone(); 
    mydiff->zero();
    for ( int i = 0; i < SampleGenerator<Real>::numMySamples(); i++ ) {
      mydiff->axpy(SampleGenerator<Real>::getMyWeight(i),(*vals[i]));
    }
    Teuchos::RCP<Vector<Real> > diff = mydiff->clone();
    SampleGenerator<Real>::sumAll(*mydiff,*diff);
    Real error = diff->norm();
    // Update global error and index sets
    error_ += error;
    if ( activeIndex_.end() != activeIndex_.begin() ) {
      activeIndex_.insert(activeIndex_.end()--,
        std::pair<Real,std::vector<int> >(error,search_index_));
    }
    else {
      activeIndex_.insert(std::pair<Real,std::vector<int> >(error,search_index_));
    }
    // Return
    return error_;
  }
  else {
    return static_cast<Real>(0);
  }
}

template<class Real>
void SparseGridGenerator<Real>::refine(void) {
  if ( adaptive_ ) {
    // Select index to investigate
    if ( !(index_.empty()) && (direction_ < (info_.dim-1)) ) {
      // Select index corresponding to next direction
      search_index_ = index_;
      direction_++;
      search_index_[direction_]++;
    }
    else {
      // Select index corresponding to largest error
      if ( index_.empty() ) {
        index_.resize(info_.dim,1);
        search_index_ = index_;
        direction_ = info_.dim;
      }
      else {
        typename std::multimap<Real,std::vector<int> >::iterator it = activeIndex_.end(); 
        it--; 
        error_ -= it->first;
        index_  = it->second;
        activeIndex_.erase(it);
        oldIndex_.insert(oldIndex_.end(),index_); 
        search_index_ = index_;
        direction_ = 0;
        search_index_[direction_]++;
      }
    }
    // Check to see if maxLevel is violated
    Quadrature<Real> rule(info_.dim);
    if ( !checkMaxLevel(search_index_) ) {
      // Check admissibility of index
      if ( isAdmissible(search_index_,direction_) ) {
        // Build difference rule
        buildDiffRule(rule,search_index_);
      }
    }
    // Set values of difference rule as points and weights
    updateSamples(rule);
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
//  if ( !(SampleGenerator<Real>::batchID()) && !inConstructor) {
//    std::cout << "Number of Global Points: " << grid_->getNumPoints() << "\n";
//  }
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

} // namespace ROL

#endif
