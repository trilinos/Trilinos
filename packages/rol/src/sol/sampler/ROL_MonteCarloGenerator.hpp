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

#ifndef ROL_MONTECARLOGENERATOR_HPP
#define ROL_MONTECARLOGENERATOR_HPP

#include "ROL_SampleGenerator.hpp"
#include "ROL_Distribution.hpp"

namespace ROL {

template<class Real>
class MonteCarloGenerator : public SampleGenerator<Real> {
private: 
  int  nSamp_;
  const bool use_normal_;
  const bool use_SA_;
  const bool adaptive_;
  const int  numNewSamps_;
  std::vector<std::vector<Real> > data_; 

  Real sum_val_;
  Real sum_val2_;
  Real sum_ng_;
  Real sum_ng2_;
  
  const bool useDist_;
  const std::vector<ROL::Ptr<Distribution<Real> > > dist_;

  Real ierf(Real input) const {
    std::vector<Real> coeff;
    Real pi = Teuchos::ScalarTraits<Real>::pi(), zero(0), one(1), two(2), tol(1e-4);
    Real c(1);
    Real tmp = c * (std::sqrt(pi)/two * input);
    Real val = tmp;
    coeff.push_back(c);
    int  cnt = 1;
    while (std::abs(tmp) > tol*std::abs(val)) {
      c = zero;
      for ( unsigned i = 0; i < coeff.size(); i++ ) {
        c += coeff[i]*coeff[coeff.size()-1-i]/((i+1)*(2*i+1));
      }
      Real ind = static_cast<Real>(cnt);
      tmp  = c/(two*ind+one) * std::pow(std::sqrt(pi)/two*input, two*ind+one);
      val += tmp;
      coeff.push_back(c);
      cnt++;
    }
    return val;
  }

  Real random(void) const {
    return static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
  }

  std::vector<std::vector<Real> > sample(int nSamp, bool store = true) {
    srand(123454321);
    const Real zero(0), one(1), two(2), tol(0.1);
    int rank = SampleGenerator<Real>::batchID();
    const int dim = (!useDist_ ? data_.size() : dist_.size());
    std::vector<Real> pts(nSamp*dim, zero);
    if (rank == 0) { 
      // Generate samples
      for (int i = 0; i < nSamp; ++i) {
        if ( !useDist_ ) {
          for (int j = 0; j < dim; ++j) {
            if ( use_normal_ ) {
              pts[j + i*dim] = std::sqrt(two*(data_[j])[1])*ierf(two*random()-one) + (data_[j])[0];
            }
            else {
              pts[j + i*dim] = ((data_[j])[1]-(data_[j])[0])*random()+(data_[j])[0];
            }
          }
        }
        else {
          for (int j = 0; j < dim; ++j) {
            pts[j + i*dim] = (dist_[j])->invertCDF(random());
            while (std::abs(pts[j + i*dim]) > tol*ROL_INF<Real>()) {
              pts[j + i*dim] = (dist_[j])->invertCDF(random());
            }
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
    for (int i = 0; i < rank; ++i) {
      offset += frac + ((i < rem) ? 1 : 0);
    }
    std::vector<std::vector<Real> > mypts;
    std::vector<Real> pt(dim);
    for (int i = 0; i < N; ++i) {
      int I = offset+i;
      for (int j = 0; j < dim; ++j) {
        pt[j] = pts[j + I*dim];
      }
      mypts.push_back(pt);
    }
    if ( store ) {
      std::vector<Real> mywts(N, one/static_cast<Real>(nSamp));
      SampleGenerator<Real>::setPoints(mypts);
      SampleGenerator<Real>::setWeights(mywts);
    }
    return mypts;
  }

  void sample(void) {
    sample(nSamp_,true);
  }

public:
  MonteCarloGenerator(const int nSamp,
                      const std::vector<ROL::Ptr<Distribution<Real> > > &dist, 
                      const ROL::Ptr<BatchManager<Real> > &bman, 
                      const bool use_SA = false,
                      const bool adaptive = false,
                      const int numNewSamps = 0)
    : SampleGenerator<Real>(bman),
      nSamp_(nSamp),
      use_normal_(false),
      use_SA_(use_SA),
      adaptive_(adaptive), 
      numNewSamps_(numNewSamps),
      sum_val_(0),
      sum_val2_(0),
      sum_ng_(0),
      sum_ng2_(0), 
      useDist_(true),
      dist_(dist) {
    int nProc = SampleGenerator<Real>::numBatches();
    TEUCHOS_TEST_FOR_EXCEPTION( nSamp_ < nProc, std::invalid_argument,
      ">>> ERROR (ROL::MonteCarloGenerator): Total number of samples is less than the number of batches!"); 
    sample();
  }

  MonteCarloGenerator(const int nSamp,
                            std::vector<std::vector<Real> > &bounds, 
                      const ROL::Ptr<BatchManager<Real> > &bman,  
                      const bool use_SA = false,
                      const bool adaptive = false,
                      const int numNewSamps = 0)
    : SampleGenerator<Real>(bman),
      nSamp_(nSamp),
      use_normal_(false),
      use_SA_(use_SA),
      adaptive_(adaptive),
      numNewSamps_(numNewSamps),
      sum_val_(0),
      sum_val2_(0),
      sum_ng_(0),
      sum_ng2_(0),
      useDist_(false) {
    int nProc = SampleGenerator<Real>::numBatches();
    TEUCHOS_TEST_FOR_EXCEPTION( nSamp_ < nProc, std::invalid_argument,
      ">>> ERROR (ROL::MonteCarloGenerator): Total number of samples is less than the number of batches!"); 
    unsigned dim = bounds.size();
    data_.clear();
    Real tmp(0);
    for ( unsigned j = 0; j < dim; j++ ) {
      if ( (bounds[j])[0] > (bounds[j])[1] ) {
        tmp = (bounds[j])[0];
        (bounds[j])[0] = (bounds[j])[1];
        (bounds[j])[1] = tmp;
        data_.push_back(bounds[j]);
      }
      data_.push_back(bounds[j]);
    }  
    sample();
  }

  MonteCarloGenerator(const int nSamp,
                      const std::vector<Real> &mean,
                      const std::vector<Real> &std, 
                      const ROL::Ptr<BatchManager<Real> > &bman,
                      const bool use_SA = false,
                      const bool adaptive = false,
                      const int numNewSamps = 0 )
    : SampleGenerator<Real>(bman),
      nSamp_(nSamp),
      use_normal_(true),
      use_SA_(use_SA),
      adaptive_(adaptive),
      numNewSamps_(numNewSamps),
      sum_val_(0),
      sum_val2_(0),
      sum_ng_(0),
      sum_ng2_(0), 
      useDist_(false) {
    int nProc = SampleGenerator<Real>::numBatches();
    TEUCHOS_TEST_FOR_EXCEPTION( nSamp_ < nProc, std::invalid_argument,
      ">>> ERROR (ROL::MonteCarloGenerator): Total number of samples is less than the number of batches!"); 
    unsigned dim = mean.size();
    data_.clear();
    std::vector<Real> tmp(2,static_cast<Real>(0));
    for ( unsigned j = 0; j < dim; j++ ) {
      tmp[0] = mean[j];
      tmp[1] = std[j];
      data_.push_back(tmp);
    }
    sample();
  }

  void update( const Vector<Real> &x ) {
    SampleGenerator<Real>::update(x);
    Real zero(0);
    sum_val_  = zero;
    sum_val2_ = zero;
    sum_ng_   = zero;
    sum_ng_   = zero;
    if ( use_SA_ ) {
      sample();
    }
  }

  Real computeError( std::vector<Real> &vals ) {
    if ( adaptive_ && !use_SA_ ) {
      Real zero(0), one(1), tol(1e-8);
      // Compute unbiased sample variance
      int cnt = 0;
      for ( int i = SampleGenerator<Real>::start(); i < SampleGenerator<Real>::numMySamples(); i++ ) {
        sum_val_  += vals[cnt];
        sum_val2_ += vals[cnt]*vals[cnt];  
        cnt++;
      }
      Real mymean = sum_val_ / nSamp_;
      Real mean   = zero;
      SampleGenerator<Real>::sumAll(&mymean,&mean,1);

      Real myvar  = (sum_val2_ - mean*mean)/(nSamp_-one);
      Real var    = zero;
      SampleGenerator<Real>::sumAll(&myvar,&var,1);
      // Return Monte Carlo error
      vals.clear();
      return std::sqrt(var/(nSamp_))*tol;
    }
    else {
      vals.clear();
      return static_cast<Real>(0);
    }
  }

  Real computeError( std::vector<ROL::Ptr<Vector<Real> > > &vals, const Vector<Real> &x ) {
    if ( adaptive_ && !use_SA_ ) {
      Real zero(0), one(1), tol(1e-4);
      // Compute unbiased sample variance
      int cnt = 0;
      Real ng = zero;
      for ( int i = SampleGenerator<Real>::start(); i < SampleGenerator<Real>::numMySamples(); i++ ) {
        ng = (vals[cnt])->norm();
        sum_ng_  += ng;
        sum_ng2_ += ng*ng;  
        cnt++;
      }
      Real mymean = sum_ng_ / nSamp_;
      Real mean   = zero;
      SampleGenerator<Real>::sumAll(&mymean,&mean,1);

      Real myvar  = (sum_ng2_ - mean*mean)/(nSamp_-one);
      Real var    = zero;
      SampleGenerator<Real>::sumAll(&myvar,&var,1);
      // Return Monte Carlo error
      vals.clear();
      return std::sqrt(var/(nSamp_))*tol;
    }
    else {
      vals.clear();
      return static_cast<Real>(0);
    }
  }

  void refine(void) {
    if ( adaptive_ && !use_SA_ ) {
      Real zero(0), one(1);
      std::vector<std::vector<Real> > pts;
      std::vector<Real> pt(data_.size(),zero);
      for ( int i = 0; i < SampleGenerator<Real>::numMySamples(); i++ ) {
        pt = SampleGenerator<Real>::getMyPoint(i);
        pts.push_back(pt);
      }
      std::vector<std::vector<Real> > pts_new = sample(numNewSamps_,false);
      pts.insert(pts.end(),pts_new.begin(),pts_new.end());
      nSamp_ += numNewSamps_;
      std::vector<Real> wts(pts.size(),one/((Real)nSamp_));
      SampleGenerator<Real>::refine();
      SampleGenerator<Real>::setPoints(pts);
      SampleGenerator<Real>::setWeights(wts);
    }
  }

};

}

#endif
