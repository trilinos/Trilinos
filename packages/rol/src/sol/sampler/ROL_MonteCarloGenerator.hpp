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
  bool use_normal_;
  bool use_SA_;
  bool adaptive_;
  int  numNewSamps_;
  std::vector<std::vector<Real> > data_; 

  Real sum_val_;
  Real sum_val2_;
  Real sum_ng_;
  Real sum_ng2_;
  
  bool useDist_;
  std::vector<Teuchos::RCP<ROL::Distribution<Real> > > dist_;

  Real ierf(Real input) {
    std::vector<Real> coeff;
    Real c   = 1.0;
    Real tmp = c * (std::sqrt(M_PI)/2.0 * input);
    Real val = tmp;
    coeff.push_back(c);
    int  cnt = 1;
    while (std::abs(tmp) > 1.e-4*std::abs(val)) {
      c = 0.0;
      for ( unsigned i = 0; i < coeff.size(); i++ ) {
        c += coeff[i]*coeff[coeff.size()-1-i]/((i+1)*(2*i+1));
      }
      tmp  = c/(2.0*(Real)cnt+1.0) * std::pow(std::sqrt(M_PI)/2.0 * input,2.0*(Real)cnt+1.0);
      val += tmp;
      coeff.push_back(c);
      cnt++;
    }
    return val;
  }

  void sample(void) {
    // Get process rank and number of processes
    int rank  = SampleGenerator<Real>::batchID();
    int nProc = SampleGenerator<Real>::numBatches();
    // Separate samples across processes
    int frac = this->nSamp_ / nProc;
    int rem  = this->nSamp_ % nProc;
    unsigned N    = frac;
    if ( rank < rem ) {
      N++;
    }
    // Generate samples
    std::vector<std::vector<Real> > pts;
    std::vector<Real> p;
    //srand((rank+1)*(rank+1)*time(NULL));
    for ( unsigned i = 0; i < N; i++ ) {
      srand(123456*(rank*N + (i+1)));
      if ( !useDist_ ) {
        p.resize(this->data_.size(),0.0);
        for ( unsigned j = 0; j < this->data_.size(); j++ ) {
          if ( this->use_normal_ ) {
            p[j] = std::sqrt(2.0*(this->data_[j])[1])*this->ierf(2.0*((Real)rand())/((Real)RAND_MAX)-1.0) + 
                   (this->data_[j])[0];
          }
          else {
            p[j] = ((this->data_[j])[1]-(this->data_[j])[0])*((Real)rand())/((Real)RAND_MAX)+(this->data_[j])[0];
          }
        }
      }
      else {
        p.resize(dist_.size(),0.0);
        for ( unsigned j = 0; j < dist_.size(); j++ ) {
          p[j] = (dist_[j])->invcdf((Real)rand()/(Real)RAND_MAX);
        }
      }
      pts.push_back(p);
    }
    std::vector<Real> wts(N,1.0/((Real)this->nSamp_));
    SampleGenerator<Real>::setPoints(pts);
    SampleGenerator<Real>::setWeights(wts);
  }

  std::vector<std::vector<Real> > sample(int nSamp, bool store = true) {
    // Get process rank and number of processes
    int rank  = SampleGenerator<Real>::batchID();
    int nProc = SampleGenerator<Real>::numBatches();
    // Separate samples across processes
    int frac = nSamp / nProc;
    int rem  = nSamp % nProc;
    unsigned N = frac;
    if ( rank < rem ) {
      N++;
    }
    // Generate samples
    std::vector<std::vector<Real> > pts;
    std::vector<Real> p;
    //srand((rank+1)*(rank+1)*time(NULL));
    for ( unsigned i = 0; i < N; i++ ) {
      srand(123456*(rank*N + (i+1)));
      if ( !useDist_ ) { 
        p.resize(this->data_.size(),0.0);
        for ( unsigned j = 0; j < this->data_.size(); j++ ) {
          if ( this->use_normal_ ) {
            p[j] = std::sqrt(2.0*(this->data_[j])[1])*this->ierf(2.0*((Real)rand())/((Real)RAND_MAX)-1.0) + 
                   (this->data_[j])[0];
          }
          else {
            p[j] = ((this->data_[j])[1]-(this->data_[j])[0])*((Real)rand())/((Real)RAND_MAX)+(this->data_[j])[0];
          }
        }
      }
      else {
        p.resize(dist_.size(),0.0);
        for ( unsigned j = 0; j < dist_.size(); j++ ) {
          p[j] = (dist_[j])->invcdf((Real)rand()/(Real)RAND_MAX);
        }
      }
      pts.push_back(p);
    }
    if ( store ) {
      std::vector<Real> wts(N,1.0/((Real)nSamp));
      SampleGenerator<Real>::setPoints(pts);
      SampleGenerator<Real>::setWeights(wts);
    }
    return pts;
  }

public:
  MonteCarloGenerator( int nSamp, std::vector<Teuchos::RCP<Distribution<Real> > > &dist, 
                       Teuchos::RCP<BatchManager<Real> > &bman, 
                       bool use_SA = false, bool adaptive = false, int numNewSamps = 0 ) 
    : SampleGenerator<Real>(bman), nSamp_(nSamp), use_normal_(false), use_SA_(use_SA), adaptive_(adaptive), 
      numNewSamps_(numNewSamps), sum_val_(0.0), sum_val2_(0.0), sum_ng_(0.0), sum_ng2_(0.0), 
      useDist_(true), dist_(dist) {
    sample();
  }

  MonteCarloGenerator( int nSamp, std::vector<std::vector<Real> > &bounds, 
                       Teuchos::RCP<BatchManager<Real> > &bman,  
                       bool use_SA = false, bool adaptive = false, int numNewSamps = 0 ) 
    : SampleGenerator<Real>(bman), nSamp_(nSamp), use_normal_(false), use_SA_(use_SA), adaptive_(adaptive),
      numNewSamps_(numNewSamps), sum_val_(0.0), sum_val2_(0.0), sum_ng_(0.0), sum_ng2_(0.0),
      useDist_(false) {
    unsigned dim = bounds.size();
    data_.clear();
    Real tmp = 0.0;
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

  MonteCarloGenerator( int nSamp, std::vector<Real> mean, std::vector<Real> std, 
                       Teuchos::RCP<BatchManager<Real> > &bman,
                       bool use_SA = false, bool adaptive = false, int numNewSamps = 0 )
    : SampleGenerator<Real>(bman), nSamp_(nSamp), use_normal_(true), use_SA_(use_SA), adaptive_(adaptive),
      numNewSamps_(numNewSamps), sum_val_(0.0), sum_val2_(0.0), sum_ng_(0.0), sum_ng2_(0.0), 
      useDist_(false) {
    unsigned dim = mean.size();
    data_.clear();
    std::vector<Real> tmp(2,0.0);
    for ( unsigned j = 0; j < dim; j++ ) {
      tmp[0] = mean[j];
      tmp[1] = std[j];
      data_.push_back(tmp);
    }
    sample();
  }

  void update( const Vector<Real> &x ) {
    SampleGenerator<Real>::update(x);
    this->sum_val_  = 0.0;
    this->sum_val2_ = 0.0;
    this->sum_ng_   = 0.0;
    this->sum_ng_   = 0.0;
    if ( this->use_SA_ ) {
      this->sample();
    }
  }

  Real computeError( std::vector<Real> &vals ) {
    if ( this->adaptive_ && !(this->use_SA_) ) {
      // Compute unbiased sample variance
      int cnt = 0;
      for ( int i = SampleGenerator<Real>::start(); i < SampleGenerator<Real>::numMySamples(); i++ ) {
        this->sum_val_  += vals[cnt];
        this->sum_val2_ += vals[cnt]*vals[cnt];  
        cnt++;
      }
      Real mymean = this->sum_val_ / this->nSamp_;
      Real mean   = 0.0;
      SampleGenerator<Real>::sumAll(&mymean,&mean,1);

      Real myvar  = (this->sum_val2_ - mean*mean)/(this->nSamp_-1.0);
      Real var    = 0.0;
      SampleGenerator<Real>::sumAll(&myvar,&var,1);
      // Return Monte Carlo error
      vals.clear();
      return std::sqrt(var/(this->nSamp_))*1.e-8;
    }
    else {
      vals.clear();
      return 0.0;
    }
  }

  Real computeError( std::vector<Teuchos::RCP<Vector<Real> > > &vals, const Vector<Real> &x ) {
    if ( this->adaptive_ && !(this->use_SA_) ) {
      // Compute unbiased sample variance
      int cnt = 0;
      Real ng = 0.0;
      for ( int i = SampleGenerator<Real>::start(); i < SampleGenerator<Real>::numMySamples(); i++ ) {
        ng = (vals[cnt])->norm();
        this->sum_ng_  += ng;
        this->sum_ng2_ += ng*ng;  
        cnt++;
      }
      Real mymean = this->sum_ng_ / this->nSamp_;
      Real mean   = 0.0;
      SampleGenerator<Real>::sumAll(&mymean,&mean,1);

      Real myvar  = (this->sum_ng2_ - mean*mean)/(this->nSamp_-1.0);
      Real var    = 0.0;
      SampleGenerator<Real>::sumAll(&myvar,&var,1);
      // Return Monte Carlo error
      vals.clear();
      return std::sqrt(var/(this->nSamp_))*1.e-4;
    }
    else {
      vals.clear();
      return 0.0;
    }
  }

  void refine(void) {
    if ( this->adaptive_ && !(this->use_SA_) ) {
      std::vector<std::vector<Real> > pts;
      std::vector<Real> pt(this->data_.size(),0.0);
      for ( int i = 0; i < SampleGenerator<Real>::numMySamples(); i++ ) {
        pt = SampleGenerator<Real>::getMyPoint(i);
        pts.push_back(pt);
      }
      std::vector<std::vector<Real> > pts_new = this->sample(this->numNewSamps_,false);
      pts.insert(pts.end(),pts_new.begin(),pts_new.end());
      this->nSamp_ += this->numNewSamps_;
      std::vector<Real> wts(pts.size(),1.0/((Real)this->nSamp_));
      SampleGenerator<Real>::refine();
      SampleGenerator<Real>::setPoints(pts);
      SampleGenerator<Real>::setWeights(wts);
    }
  }

};

}

#endif
