// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_MONTECARLOGENERATORDEF_HPP
#define ROL_MONTECARLOGENERATORDEF_HPP

namespace ROL {

template<typename Real>
Real MonteCarloGenerator<Real>::ierf(Real input) const {
  std::vector<Real> coeff;
  const Real pi = ScalarTraits<Real>::pi(), zero(0), one(1), two(2), tol(1e-4);
  Real c(1);
  Real tmp = c * (std::sqrt(pi)/two * input);
  Real val = tmp;
  coeff.push_back(c);
  int  cnt = 1;
  while (std::abs(tmp) > tol*std::abs(val)) {
    c = zero;
    for ( unsigned i = 0; i < coeff.size(); ++i )
      c += coeff[i]*coeff[coeff.size()-1-i]/static_cast<Real>((i+1)*(2*i+1));
    Real ind = static_cast<Real>(cnt);
    tmp  = c/(two*ind+one) * std::pow(std::sqrt(pi)/two*input, two*ind+one);
    val += tmp;
    coeff.push_back(c);
    cnt++;
  }
  return val;
}

template<typename Real>
Real MonteCarloGenerator<Real>::random(void) const {
  return static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
}

template<typename Real>
std::vector<std::vector<Real>> MonteCarloGenerator<Real>::sample(int nSamp, bool store, bool refine) {
  if (!refine) srand(seed_);
  const Real zero(0), one(1), two(2), tol(0.1);
  int rank = SampleGenerator<Real>::batchID();
  const int dim = (!useDist_ ? data_.size() : dist_.size());
  std::vector<Real> pts(nSamp*dim, zero);
  if (rank == 0) { 
    // Generate samples
    for (int i = 0; i < nSamp; ++i) {
      if ( !useDist_ ) {
        for (int j = 0; j < dim; ++j) {
          if ( use_normal_ )
            pts[j + i*dim] = std::sqrt(two*(data_[j])[1])*ierf(two*random()-one) + (data_[j])[0];
          else
            pts[j + i*dim] = ((data_[j])[1]-(data_[j])[0])*random()+(data_[j])[0];
        }
      }
      else {
        for (int j = 0; j < dim; ++j) {
          pts[j + i*dim] = (dist_[j])->invertCDF(random());
          while (std::abs(pts[j + i*dim]) > tol*ROL_INF<Real>())
            pts[j + i*dim] = (dist_[j])->invertCDF(random());
        }
      }
    }
  }
  SampleGenerator<Real>::broadcast(&pts[0],nSamp*dim,0);
  // Separate samples across processes
  int nProc  = SampleGenerator<Real>::numBatches();
  int frac   = nSamp / nProc;
  int rem    = nSamp % nProc;
  int N      = frac + ((rank < rem) ? 1 : 0);
  int offset = 0;
  for (int i = 0; i < rank; ++i) offset += frac + ((i < rem) ? 1 : 0);
  std::vector<std::vector<Real>> mypts;
  std::vector<Real> pt(dim);
  for (int i = 0; i < N; ++i) {
    int I = offset+i;
    for (int j = 0; j < dim; ++j) pt[j] = pts[j + I*dim];
    mypts.push_back(pt);
  }
  if ( store ) {
    std::vector<Real> mywts(N, one/static_cast<Real>(nSamp));
    SampleGenerator<Real>::setPoints(mypts);
    SampleGenerator<Real>::setWeights(mywts);
  }
  return mypts;
}

template<typename Real>
MonteCarloGenerator<Real>::MonteCarloGenerator(int nSamp,
                    const std::vector<Ptr<Distribution<Real>>> &dist, 
                    const Ptr<BatchManager<Real>> &bman, 
                    bool use_SA, bool adaptive, int numNewSamps, int seed)
  : SampleGenerator<Real>(bman),
    dist_(dist),
    nSamp_(nSamp),
    use_normal_(false),
    use_SA_(use_SA),
    adaptive_(adaptive),
    numNewSamps_(numNewSamps),
    useDist_(true),
    seed_(seed),
    sum_val_(0),
    sum_val2_(0),
    sum_ng_(0),
    sum_ng2_(0) {
  int nProc = SampleGenerator<Real>::numBatches();
  ROL_TEST_FOR_EXCEPTION( nSamp_ < nProc, std::invalid_argument,
    ">>> ERROR (ROL::MonteCarloGenerator): Total number of samples is less than the number of batches!");
  sample(nSamp_,true,false);
}

template<typename Real>
MonteCarloGenerator<Real>::MonteCarloGenerator(int nSamp,
                    std::vector<std::vector<Real>> &bounds,
                    const Ptr<BatchManager<Real>> &bman,
                    bool use_SA, bool adaptive, int numNewSamps, int seed)
  : SampleGenerator<Real>(bman),
    nSamp_(nSamp),
    use_normal_(false),
    use_SA_(use_SA),
    adaptive_(adaptive),
    numNewSamps_(numNewSamps),
    useDist_(false),
    seed_(seed),
    sum_val_(0),
    sum_val2_(0),
    sum_ng_(0),
    sum_ng2_(0) {
  int nProc = SampleGenerator<Real>::numBatches();
  ROL_TEST_FOR_EXCEPTION( nSamp_ < nProc, std::invalid_argument,
    ">>> ERROR (ROL::MonteCarloGenerator): Total number of samples is less than the number of batches!");
  unsigned dim = bounds.size();
  data_.clear();
  Real tmp(0);
  for ( unsigned j = 0; j < dim; ++j ) {
    if ( (bounds[j])[0] > (bounds[j])[1] ) {
      tmp = (bounds[j])[0];
      (bounds[j])[0] = (bounds[j])[1];
      (bounds[j])[1] = tmp;
      data_.push_back(bounds[j]);
    }
    data_.push_back(bounds[j]);
  }
  sample(nSamp_,true,false);
}

template<typename Real>
MonteCarloGenerator<Real>::MonteCarloGenerator(int nSamp,
                    const std::vector<Real> &mean,
                    const std::vector<Real> &std,
                    const Ptr<BatchManager<Real>> &bman,
                    bool use_SA, bool adaptive, int numNewSamps, int seed)
  : SampleGenerator<Real>(bman),
    nSamp_(nSamp),
    use_normal_(true),
    use_SA_(use_SA),
    adaptive_(adaptive),
    numNewSamps_(numNewSamps),
    useDist_(false),
    seed_(seed),
    sum_val_(0),
    sum_val2_(0),
    sum_ng_(0),
    sum_ng2_(0) {
  int nProc = SampleGenerator<Real>::numBatches();
  ROL_TEST_FOR_EXCEPTION( nSamp_ < nProc, std::invalid_argument,
    ">>> ERROR (ROL::MonteCarloGenerator): Total number of samples is less than the number of batches!");
  unsigned dim = mean.size();
  data_.clear();
  for ( unsigned j = 0; j < dim; ++j )
    data_.push_back({mean[j],std[j]});
  sample(nSamp_,true,false);
}

template<typename Real>
void MonteCarloGenerator<Real>::update( const Vector<Real> &x ) {
  SampleGenerator<Real>::update(x);
  sum_val_  = static_cast<Real>(0);
  sum_val2_ = static_cast<Real>(0);
  sum_ng_   = static_cast<Real>(0);
  sum_ng_   = static_cast<Real>(0);
  if ( use_SA_ ) sample(nSamp_,true,true);
}

template<typename Real>
Real MonteCarloGenerator<Real>::computeError( std::vector<Real> &vals ) {
  Real err(0);
  if ( adaptive_ && !use_SA_ ) {
    const Real zero(0), one(1), tol(1e-8);
    // Compute unbiased sample variance
    int cnt = 0;
    for ( int i = SampleGenerator<Real>::start(); i < SampleGenerator<Real>::numMySamples(); ++i ) {
      sum_val_  += vals[cnt];
      sum_val2_ += vals[cnt]*vals[cnt];
      cnt++;
    }
    Real mymean = sum_val_ / static_cast<Real>(nSamp_);
    Real mean   = zero;
    SampleGenerator<Real>::sumAll(&mymean,&mean,1);

    Real myvar  = (sum_val2_ - mean*mean)/(static_cast<Real>(nSamp_)-one);
    Real var    = zero;
    SampleGenerator<Real>::sumAll(&myvar,&var,1);
    // Return Monte Carlo error
    vals.clear();
    err = std::sqrt(var/static_cast<Real>(nSamp_))*tol;
  }
  else {
    vals.clear();
  }
  return err;
}

template<typename Real>
Real MonteCarloGenerator<Real>::computeError( std::vector<Ptr<Vector<Real>>> &vals,
                                              const Vector<Real> &x ) {
  Real err(0);
  if ( adaptive_ && !use_SA_ ) {
    const Real zero(0), one(1), tol(1e-4);
    // Compute unbiased sample variance
    int cnt = 0;
    Real ng = zero;
    for ( int i = SampleGenerator<Real>::start(); i < SampleGenerator<Real>::numMySamples(); ++i ) {
      ng = (vals[cnt])->norm();
      sum_ng_  += ng;
      sum_ng2_ += ng*ng;
      cnt++;
    }
    Real mymean = sum_ng_ / static_cast<Real>(nSamp_);
    Real mean   = zero;
    SampleGenerator<Real>::sumAll(&mymean,&mean,1);

    Real myvar  = (sum_ng2_ - mean*mean)/(static_cast<Real>(nSamp_)-one);
    Real var    = zero;
    SampleGenerator<Real>::sumAll(&myvar,&var,1);
    // Return Monte Carlo error
    vals.clear();
    err = std::sqrt(var/static_cast<Real>(nSamp_))*tol;
  }
  else {
    vals.clear();
  }
  return err;
}

template<typename Real>
void MonteCarloGenerator<Real>::refine(void) {
  if ( adaptive_ && !use_SA_ ) {
    std::vector<std::vector<Real>> pts;
    for ( int i = 0; i < SampleGenerator<Real>::numMySamples(); ++i )
      pts.push_back(SampleGenerator<Real>::getMyPoint(i));
    std::vector<std::vector<Real>> pts_new = sample(numNewSamps_,false,true);
    pts.insert(pts.end(),pts_new.begin(),pts_new.end());
    nSamp_ += numNewSamps_;
    std::vector<Real> wts(pts.size(),static_cast<Real>(1)/static_cast<Real>(nSamp_));
    SampleGenerator<Real>::refine();
    SampleGenerator<Real>::setPoints(pts);
    SampleGenerator<Real>::setWeights(wts);
  }
}

template<typename Real>
int MonteCarloGenerator<Real>::numGlobalSamples(void) const {
  return nSamp_;
}

}
#endif
