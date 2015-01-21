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
SparseGridGenerator<Real>::SparseGridGenerator( Teuchos::RCP<BatchManager<Real> > &bman, 
                                                SparseGridInfo &info, bool adaptive )
  : SampleGenerator<Real>(bman), adaptive_(adaptive), info_(info) {
  if ( adaptive ) {
    //index_.resize(info.dim,1);
    index_.clear();
    direction_ = 0;
    error_ = 0.0;
    activeIndex_.clear();
    oldIndex_.clear();
    grid_ = Teuchos::rcp(new Quadrature<Real>(info.dim));
  }
  else {
    int dim_num   = info.dim; 
    int level_max = info.maxLevel - 1;
    std::vector<EROLBurkardt> rule1D = info.rule1D;
    std::vector<EROLGrowth> growth1D = info.growth1D;
    bool isNormalized = info.normalize;
    bool useSandia = info.useSandia;
    grid_ = Teuchos::rcp(new Quadrature<Real>(dim_num,level_max,rule1D,growth1D,isNormalized,useSandia));
    // Split points and weights across processors
    this->setSamples(true);
  }
}

template <class Real>
SparseGridGenerator<Real>::SparseGridGenerator(Teuchos::RCP<BatchManager<Real> > &bman, 
                                               const char* SGinfo,const char* SGdata, 
                                               bool isNormalized) : SampleGenerator<Real>(bman), 
                                                                    adaptive_(false) {
  grid_ = Teuchos::rcp(new Quadrature<Real>(SGinfo,SGdata,isNormalized));
  // Split points and weights across processors
  this->setSamples(true);
}

template<class Real>
bool SparseGridGenerator<Real>::checkMaxLevel(std::vector<int> &index) {
  int  level = 0;
  bool useMax = true;
  for ( int l = 0; l < index.size(); l++ ) {
    if ( useMax ) {
      level = std::max(level,index[l]);
    }
    else {
      level += index[l];
    }
  }
  if ( level >= (this->info_).maxLevel ) {
    return true;
  }
  return false;
}

template<class Real>
bool SparseGridGenerator<Real>::isAdmissible(std::vector<int> &index, int direction) {
  for ( int i = 0; i < (this->info_).dim; i++ ) {
    if ( index[i] > 1 && i != direction ) {
      index[i]--;
      if ( !((this->oldIndex_).count(index)) ) {
        return false;
      }
      index[i]++;
    }
  }
  return true;
}

template<class Real> 
void SparseGridGenerator<Real>::buildDiffRule(Quadrature<Real> &outRule, std::vector<int> &index) {
  int numPoints  = 0;
  for ( int i = 0; i < (this->info_).dim; i++ ) {
    numPoints = growthRule1D(index[i],((this->info_).growth1D)[i],((this->info_).rule1D)[i]);
    Quadrature<Real> diffRule(((this->info_).rule1D)[i],numPoints,((this->info_).normalize));
    if ( numPoints != 1 ) {
      numPoints = growthRule1D(index[i]-1,((this->info_).growth1D)[i],((this->info_).rule1D)[i]);
      Quadrature<Real> rule(((this->info_).rule1D)[i],numPoints,((this->info_).normalize));
      diffRule.update(-1.0,rule);
    }
    outRule = kron_prod<Real>(outRule,diffRule);
  }
}

template<class Real>
void SparseGridGenerator<Real>::update(const Vector<Real> &x) {
  SampleGenerator<Real>::update(x);
  if ( this->adaptive_ ) {
    //(this->index_).resize((this->info_).dim,1);
    (this->index_).clear();
    this->direction_ = 0;
    this->error_ = 0.0;
    (this->activeIndex_).clear();
    (this->oldIndex_).clear();
    this->grid_ = Teuchos::rcp(new Quadrature<Real>((this->info_).dim));
  }
}

template<class Real>
Real SparseGridGenerator<Real>::computeError(std::vector<Real> &vals){
  if ( this->adaptive_ ) {
    Real myerror = 0.0;
    for ( int i = 0; i < SampleGenerator<Real>::numMySamples(); i++ ) {
      myerror += vals[i]*SampleGenerator<Real>::getMyWeight(i);
    }
    Real error = 0.0;
    SampleGenerator<Real>::sumAll(&myerror,&error,1);
    error = std::abs(error);
    // Update global error and index sets
    this->error_ += error;
    if ( (this->activeIndex_).end() != (this->activeIndex_).begin() ) {
      (this->activeIndex_).insert((this->activeIndex_).end()--,
        std::pair<Real,std::vector<int> >(error,this->search_index_));
    }
    else {
      (this->activeIndex_).insert(std::pair<Real,std::vector<int> >(error,this->search_index_));
    }
    // Return
    return this->error_;
  }
  else {
    return 0.0;
  }
}

template<class Real>
Real SparseGridGenerator<Real>::computeError(std::vector<Teuchos::RCP<Vector<Real> > > &vals, const Vector<Real> &x ){
  if ( this->adaptive_ ) {
    Teuchos::RCP<Vector<Real> > mydiff = x.clone(); 
    mydiff->zero();
    for ( int i = 0; i < SampleGenerator<Real>::numMySamples(); i++ ) {
      mydiff->axpy(SampleGenerator<Real>::getMyWeight(i),(*vals[i]));
    }
    Teuchos::RCP<Vector<Real> > diff = mydiff->clone();
    SampleGenerator<Real>::sumAll(*mydiff,*diff);
    Real error = diff->norm();
    // Update global error and index sets
    this->error_ += error;
    if ( (this->activeIndex_).end() != (this->activeIndex_).begin() ) {
      (this->activeIndex_).insert((this->activeIndex_).end()--,
        std::pair<Real,std::vector<int> >(error,this->search_index_));
    }
    else {
      (this->activeIndex_).insert(std::pair<Real,std::vector<int> >(error,this->search_index_));
    }
    // Return
    return this->error_;
  }
  else {
    return 0.0;
  }
}

template<class Real>
void SparseGridGenerator<Real>::refine(void) {
  if ( this->adaptive_ ) {
    // Select index to investigate
    if ( !((this->index_).empty()) && ((this->direction_) < ((this->info_).dim-1)) ) {
      // Select index corresponding to next direction
      this->search_index_ = (this->index_);
      this->direction_++;
      (this->search_index_)[this->direction_]++;
    }
    else {
      // Select index corresponding to largest error
      if ( (this->index_).empty() ) {
        (this->index_).resize((this->info_).dim,1);
        this->search_index_ = this->index_;
        this->direction_ = (this->info_).dim;
      }
      else {
        typename std::multimap<Real,std::vector<int> >::iterator it = (this->activeIndex_).end(); 
        it--; 
        this->error_ -= it->first;
        this->index_  = it->second;
        (this->activeIndex_).erase(it);
        (this->oldIndex_).insert((this->oldIndex_).end(),(this->index_)); 
        this->search_index_ = (this->index_);
        this->direction_ = 0;
        (this->search_index_)[this->direction_]++;
      }
    }
    // Check to see if maxLevel is violated
    Quadrature<Real> rule((this->info_).dim);
    if ( !(this->checkMaxLevel(this->search_index_)) ) {
      // Check admissibility of index
      if ( this->isAdmissible((this->search_index_),this->direction_) ) {
        // Build difference rule
        this->buildDiffRule(rule,this->search_index_);
      }
    }
    // Set values of difference rule as points and weights
    this->updateSamples(rule);
  }
}

template<class Real> 
void SparseGridGenerator<Real>::splitSamples(std::vector<std::vector<Real> > &mypts, std::vector<Real> &mywts) {
  // Get global points and weights
  std::vector<std::vector<Real> > pts;
  std::vector<Real> wts;
  this->grid_->getCubature(pts,wts);
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
  if ( !(SampleGenerator<Real>::batchID()) && !inConstructor) {
    std::cout << "Number of Global Points: " << (this->grid_)->getNumPoints() << "\n";
  }
  if ( this->adaptive_ || inConstructor ) {
    // Split samples based on PID
    std::vector<std::vector<Real> > mypts;
    std::vector<Real> mywts;
    this->splitSamples(mypts,mywts);
    // Set local points and weights
    SampleGenerator<Real>::setPoints(mypts);
    SampleGenerator<Real>::setWeights(mywts);
  }
}

template<class Real>
void SparseGridGenerator<Real>::updateSamples( Quadrature<Real> &grid ) {
  // Add increment to stored sparse grid
  this->grid_->update(1.0,grid);
  // Split global stored points and weights across processors
  std::vector<std::vector<Real> > mypts;
  std::vector<Real> mywts;
  this->splitSamples(mypts,mywts);
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
