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

#ifndef ROL_SROMGENERATOR_HPP
#define ROL_SROMGENERATOR_HPP

#include "ROL_SampleGenerator.hpp"

#include "ROL_Objective.hpp"
#include "ROL_BoundConstraint.hpp"
#include "ROL_EqualityConstraint.hpp"

#include "ROL_Algorithm.hpp"
#include "ROL_BoundConstraint.hpp"

#include "ROL_MomentObjective.hpp"
#include "ROL_CDFObjective.hpp"
#include "ROL_LinearCombinationObjective.hpp"
//#include "ROL_SROMBoundConstraint.hpp"
#include "ROL_SROMEqualityConstraint.hpp"
#include "ROL_SROMVector.hpp"

#include "ROL_StdVector.hpp"

namespace ROL {

template<class Real>
class SROMGenerator : public SampleGenerator<Real> {
private:
  // Parameterlist for optimization
  Teuchos::ParameterList parlist_;
  // Vector of distributions (size = dimension of space)
  std::vector<Teuchos::RCP<Distribution<Real> > > dist_;
  // Objective function (linear combination of multiple)
  Teuchos::RCP<Objective<Real> > obj_;
  std::vector<Teuchos::RCP<Objective<Real> > > obj_vec_;
  // Bound constraint information
  Teuchos::RCP<BoundConstraint<Real> > bnd_;
  Teuchos::RCP<Vector<Real> > bnd_lo_;
  Teuchos::RCP<Vector<Real> > bnd_hi_;
  // Sum of probabilities must be one
  Teuchos::RCP<EqualityConstraint<Real> > con_;
  // Optimization algorithm
  Teuchos::RCP<Algorithm<Real> > algo_;

  const size_t dimension_;
  size_t nSamp_;
  size_t numNewSamps_;
  bool adaptive_;

  void buildOptimizer(const bool useAugLag = true) {
    if ( !useAugLag ) {
      algo_ = Teuchos::rcp(new Algorithm<Real>("Moreau-Yosida Penalty",parlist_,false));
    }
    else {
      algo_ = Teuchos::rcp(new Algorithm<Real>("Augmented Lagrangian",parlist_,false));
    }
  }

  void initialize(Teuchos::RCP<Vector<Real> >       &x,
                  Teuchos::RCP<BatchManager<Real> > &bman) {
    std::vector<Real> pt(dimension_*nSamp_,0.), wt(nSamp_,1./(Real)nSamp_);
    std::vector<Real> typx(dimension_*nSamp_,1.), typw(nSamp_,1.);
    std::vector<Real> pt_lo(dimension_*nSamp_,0.), pt_hi(dimension_*nSamp_,0.);
    std::vector<Real> wt_lo(nSamp_,0.), wt_hi(nSamp_,1.);
    Real mean = 1., lo = 0., hi = 0.;
    for ( size_t j = 0; j < dimension_; j++) {
      mean = std::abs(dist_[j]->moment(1));
      mean = ((mean > ROL_EPSILON) ? mean : 1.);
      lo   = dist_[j]->lowerBound();
      hi   = dist_[j]->upperBound();
      for (size_t i = 0; i < nSamp_; i++) {
        pt[i*dimension_ + j] = dist_[j]->invertCDF((Real)rand()/(Real)RAND_MAX);
        typx[i*dimension_ + j] = mean;
        pt_lo[i*dimension_ + j] = lo;
        pt_hi[i*dimension_ + j] = hi;
      }
    }
    x = Teuchos::rcp(new PrimalSROMVector<Real>(Teuchos::rcp(new std::vector<Real>(pt)),
                                                Teuchos::rcp(new std::vector<Real>(wt)),
                                                bman,typx,typw));
    bnd_lo_ = Teuchos::rcp(new PrimalSROMVector<Real>(Teuchos::rcp(new std::vector<Real>(pt_lo)),
                                                      Teuchos::rcp(new std::vector<Real>(wt_lo)),
                                                      bman,typx,typw));
    bnd_hi_ = Teuchos::rcp(new PrimalSROMVector<Real>(Teuchos::rcp(new std::vector<Real>(pt_hi)),
                                                      Teuchos::rcp(new std::vector<Real>(wt_hi)),
                                                      bman,typx,typw));
    bnd_ = Teuchos::rcp(new BoundConstraint<Real>(bnd_lo_,bnd_hi_));
    con_ = Teuchos::rcp(new SROMEqualityConstraint<Real>(bman));
  }

  void pruneSamples(const Vector<Real> &x) {
    const PrimalSROMVector<Real> &ex = Teuchos::dyn_cast<const PrimalSROMVector<Real> >(x);
    // Remove points with zero weight
    std::vector<std::vector<Real> > pts;
    std::vector<Real> wts;
    for (size_t i = 0; i < nSamp_; i++) {
      if ( ex.getWeight(i) > ROL_EPSILON ) {
        pts.push_back(*(ex.getPoint(i)));
        wts.push_back(ex.getWeight(i));
      }
    }
    nSamp_ = wts.size();
    SampleGenerator<Real>::setPoints(pts);
    SampleGenerator<Real>::setWeights(wts);
  }

public:

  SROMGenerator(Teuchos::ParameterList                                &parlist,
                Teuchos::RCP<BatchManager<Real> >                     &bman,
                const std::vector<Teuchos::RCP<Distribution<Real> > > &dist)
    : SampleGenerator<Real>(bman), parlist_(parlist), dist_(dist),
      dimension_(dist.size()) {
    // Get SROM sublist
    Teuchos::ParameterList list = parlist.sublist("SOL").sublist("Sample Generator").sublist("SROM");
    obj_vec_.clear(); obj_vec_.resize(2);
    size_t nSamp = list.get("Number of Samples",50);
    adaptive_    = list.get("Adaptive Sampling",false);
    numNewSamps_ = list.get("Number of New Samples Per Adaptiation",0);
    // Compute batch local number of samples
    size_t rank  = (size_t)SampleGenerator<Real>::batchID();
    size_t nProc = (size_t)SampleGenerator<Real>::numBatches();
    size_t frac  = nSamp / nProc;
    size_t rem   = nSamp % nProc;
    nSamp_       = frac + ((rank < rem) ? 1 : 0);
    // Build CDF objective function
    Real scale = list.get("CDF Smoothing Parameter",1.e-2);
    std::vector<Real> pt_lo(dimension_,0.), pt_hi(dimension_,0.);
    for ( size_t j = 0; j < dimension_; j++) {
      pt_lo[j] = dist[j]->lowerBound();
      pt_hi[j] = dist[j]->upperBound();
    }
    obj_vec_[0] = Teuchos::rcp(new CDFObjective<Real>(dist,pt_lo,pt_hi,bman,scale)); 
    // Build moment matching objective function
    Teuchos::Array<int> tmp_order
      = Teuchos::getArrayFromStringParameter<int>(list,"Moments");
    std::vector<size_t> order(tmp_order.size(),0);
    for ( int i = 0; i < tmp_order.size(); i++) {
      order[i] = static_cast<size_t>(tmp_order[i]);
    }
    obj_vec_[1] = Teuchos::rcp(new MomentObjective<Real>(dist,order,bman));
    // Build linear combination objective function
    Teuchos::Array<Real> tmp_coeff
      = Teuchos::getArrayFromStringParameter<Real>(list,"Coefficients");
    std::vector<Real> coeff(2,0.);
    coeff[0] = tmp_coeff[0]; coeff[1] = tmp_coeff[1];
    obj_ = Teuchos::rcp(new LinearCombinationObjective<Real>(coeff,obj_vec_));
    // Initialize constraints and initial guess
    Teuchos::RCP<Vector<Real> > x;
    initialize(x,bman);
    StdVector<Real> l(Teuchos::rcp(new std::vector<Real>(1,0.)));
    // Solve optimization problems to sample
    bool useAugLag = true;
    buildOptimizer(useAugLag);
    algo_->run(*x,l,*obj_,*con_,*bnd_,!SampleGenerator<Real>::batchID());
    // Prune samples with zero weight and set samples/weights
    pruneSamples(*x);
  }

  void refine(void) {}
};

}

#endif

//  SROMGenerator(Teuchos::ParameterList &parlist,
//                Teuchos::RCP<BatchManager<Real> > &bman, 
//                Teuchos::RCP<Objective<Real> > &obj,
//                const size_t dimension,
//                const size_t nSamp = 10,
//                const bool adaptive = false,
//                const size_t numNewSamps = 0 ) 
//    : SampleGenerator<Real>(bman), parlist_(parlist),
//      obj_(obj), dimension_(dimension),
//      nSamp_(nSamp), numNewSamps_(numNewSamps), adaptive_(adaptive) {
//    // Build ROL algorithm and solve SROM optimization problem
//    SROMVector<Real> x(Teuchos::rcp(new std::vector<Real>(dimension_*nSamp_,0.)),
//                       Teuchos::rcp(new std::vector<Real>(nSamp_,0.)));
//    StdVector<Real> l(Teuchos::rcp(new std::vector<Real>(1,0.)));
//    //bnd_ = Teuchos::rcp(new SROMBoundConstraint<Real>(dimension_));
//    bnd_ = Teuchos::rcp(new BoundConstraint<Real>);
//    con_ = Teuchos::rcp(new SROMEqualityConstraint<Real>);
//    bool useAugLag = true;
//    buildOptimizer(useAugLag);
//    algo_->run(x,l,*obj_,*con_,*bnd_,!SampleGenerator<Real>::batchID());
//    // Prune samples with zero weight and set samples/weights
//    std::vector<std::vector<Real> > allPoints;
//    std::vector<Real> allWeights;
//    pruneSamples(allPoints,allWeights,x);
//    splitSamples(allPoints,allWeights);
//  }
//
//  SROMGenerator(Teuchos::ParameterList &parlist,
//                Teuchos::RCP<BatchManager<Real> > &bman, 
//                Teuchos::RCP<Objective<Real> > &obj,
//                Teuchos::RCP<BoundConstraint<Real> > &bnd,
//                const size_t dimension,
//                const size_t nSamp = 10,
//                const bool adaptive = false,
//                const size_t numNewSamps = 0 ) 
//    : SampleGenerator<Real>(bman), parlist_(parlist),
//      obj_(obj), bnd_(bnd), dimension_(dimension),
//      nSamp_(nSamp), numNewSamps_(numNewSamps), adaptive_(adaptive) {
//    // Build ROL algorithm and solve SROM optimization problem
//    SROMVector<Real> x(Teuchos::rcp(new std::vector<Real>(dimension_*nSamp_,0.)),
//                       Teuchos::rcp(new std::vector<Real>(nSamp_,0.)));
//    StdVector<Real> l(Teuchos::rcp(new std::vector<Real>(1,0.)));
//    con_ = Teuchos::rcp(new SROMEqualityConstraint<Real>);
//    bool useAugLag = true;
//    buildOptimizer(useAugLag);
//    algo_->run(x,l,*obj_,*con_,*bnd_,!SampleGenerator<Real>::batchID());
//    // Prune samples with zero weight and set samples/weights
//    std::vector<std::vector<Real> > allPoints;
//    std::vector<Real> allWeights;
//    pruneSamples(allPoints,allWeights,x);
//    splitSamples(allPoints,allWeights);
//  }
//
//  SROMGenerator(Teuchos::ParameterList &parlist,
//                Teuchos::RCP<BatchManager<Real> > &bman, 
//                Teuchos::RCP<Objective<Real> > &obj,
//                Teuchos::RCP<BoundConstraint<Real> > &bnd,
//                Teuchos::RCP<Vector<Real> > &x,
//                const size_t dimension,
//                const size_t nSamp = 10,
//                const bool adaptive = false,
//                const size_t numNewSamps = 0 ) 
//    : SampleGenerator<Real>(bman), parlist_(parlist), obj_(obj), bnd_(bnd),
//      dimension_(dimension), nSamp_(nSamp), numNewSamps_(numNewSamps),
//      adaptive_(adaptive) {
//    // Build ROL algorithm and solve SROM optimization problem
//    StdVector<Real> l(Teuchos::rcp(new std::vector<Real>(1,0.)));
//    con_ = Teuchos::rcp(new SROMEqualityConstraint<Real>);
//    bool useAugLag = true;
//    buildOptimizer(useAugLag);
//    algo_->run(*x,l,*obj_,*con_,*bnd_,!SampleGenerator<Real>::batchID());
//    // Prune samples with zero weight and set samples/weights
//    const SROMVector<Real> &ex = Teuchos::static_cast<const SROMVector<Real> >(*x);
//    std::vector<std::vector<Real> > allPoints;
//    std::vector<Real> allWeights;
//    pruneSamples(allPoints,allWeights,ex);
//    splitSamples(allPoints,allWeights);
//  }
