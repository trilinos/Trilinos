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

#include "ROL_MomentObjective.hpp"
#include "ROL_SROMBoundConstraint.hpp"
#include "ROL_SROMEqualityConstraint.hpp"
#include "ROL_SROMVector.hpp"

#include "ROL_StdVector.hpp"

namespace ROL {

template<class Real>
class SROMGenerator : public SampleGenerator<Real> {
private: 
  Teuchos::RCP<Objective<Real> > obj_;
  Teuchos::RCP<BoundConstraint<Real> > bnd_;
  Teuchos::RCP<EqualityConstraint<Real> > con_;

  Teuchos::RCP<DefaultAlgorithm<Real> > algo_;

  Teuchos::ParameterList parlist_;

  const size_t dimension_;
  size_t nSamp_;
  const size_t numNewSamps_;
  const bool adaptive_;

  void buildOptimizer(const bool useAugLag = true) {
    if ( !useAugLag ) {
      parlist_.sublist("Step").sublist("Moreau-Yosida Penalty").set("Initial Penalty Parameter",10.);
      parlist_.sublist("Step").sublist("Moreau-Yosida Penalty").set("Penalty Parameter Growth Factor",1.);
      parlist_.sublist("Step").sublist("Moreau-Yosida Penalty").sublist("Subproblem").set("Optimality Tolerance",1.e-8);
      parlist_.sublist("Step").sublist("Moreau-Yosida Penalty").sublist("Subproblem").set("Feasibility Tolerance",1.e-8);
      parlist_.sublist("Step").sublist("Moreau-Yosida Penalty").sublist("Subproblem").set("Iteration Limit",1000);
      parlist_.sublist("Step").sublist("Moreau-Yosida Penalty").sublist("Subproblem").set("Print History",false);

      parlist_.sublist("Status Test").get("Gradient Tolerance",   1.e-8);
      parlist_.sublist("Status Test").get("Constraint Tolerance", 1.e-8);

      algo_ = Teuchos::rcp(new DefaultAlgorithm<Real>("Moreau-Yosida Penalty",parlist_,false));
    }
    else {
      parlist_.sublist("Step").sublist("Augmented Lagrangian").get("Initial Penalty Parameter",1.e1);
      parlist_.sublist("Step").sublist("Augmented Lagrangian").get("Penalty Parameter Growth Factor",1.e2);
      parlist_.sublist("Step").sublist("Augmented Lagrangian").get("Optimality Tolerance Update Exponent",1.);
      parlist_.sublist("Step").sublist("Augmented Lagrangian").get("Feasibility Tolerance Update Exponent",0.1);
      parlist_.sublist("Step").sublist("Augmented Lagrangian").get("Optimality Tolerance Decrease Exponent",1.);
      parlist_.sublist("Step").sublist("Augmented Lagrangian").get("Feasibility Tolerance Decrease Exponent",0.9);
      parlist_.sublist("Step").sublist("Augmented Lagrangian").get("Initial Optimality Tolerance",1.);
      parlist_.sublist("Step").sublist("Augmented Lagrangian").get("Initial Feasibility Tolerance",1.);
      parlist_.sublist("Step").sublist("Augmented Lagrangian").get("Minimum Penalty Parameter Reciprocal",0.1);
      parlist_.sublist("Step").sublist("Augmented Lagrangian").get("Print Intermediate Optimization History",false);
      parlist_.sublist("Step").sublist("Augmented Lagrangian").get("Subproblem Iteration Limit",1000);
      parlist_.sublist("Step").sublist("Augmented Lagrangian").get("Subproblem Step Type",0);

      parlist_.sublist("Status Test").get("Gradient Tolerance",   1.e-8);
      parlist_.sublist("Status Test").get("Constraint Tolerance", 1.e-8);

      algo_ = Teuchos::rcp(new DefaultAlgorithm<Real>("Augmented Lagrangian",parlist_,false));
    }
  }

  void splitSamples(const std::vector<std::vector<Real> > &allPoints,
                    const std::vector<Real> &allWeights) {
    // Separate samples/weights across batches
    size_t rank  = (size_t)SampleGenerator<Real>::batchID();
    size_t nProc = (size_t)SampleGenerator<Real>::numBatches();
    size_t frac  = nSamp_ / nProc;
    size_t rem   = nSamp_ % nProc;
    size_t N     = frac + ((rank < rem) ? 1 : 0);
    size_t index = 0;

    std::vector<std::vector<Real> > pts;
    std::vector<Real> wts;
    for (size_t i = 0; i < N; i++) {
      index = i*nProc + rank;
      pts.push_back(allPoints[index]);
      wts.push_back(allWeights[index]);
    }
    SampleGenerator<Real>::setPoints(pts);
    SampleGenerator<Real>::setWeights(wts);
  }

  void pruneSamples(std::vector<std::vector<Real> > &pts, std::vector<Real> &wts,
              const SROMVector<Real> &x) {
    // Remove points with zero weight
    for (size_t i = 0; i < nSamp_; i++) {
      if ( x.getWeight(i) > ROL_EPSILON ) {
        pts.push_back(*(x.getPoint(i)));
        wts.push_back(x.getWeight(i));
      }
    }
    nSamp_ = wts.size();
  }

public:
  SROMGenerator(Teuchos::RCP<BatchManager<Real> > &bman, 
                Teuchos::RCP<Objective<Real> > &obj,
                const size_t dimension,
                const size_t nSamp = 10,
                const bool adaptive = false,
                const size_t numNewSamps = 0 ) 
    : SampleGenerator<Real>(bman), obj_(obj), dimension_(dimension),
      nSamp_(nSamp), numNewSamps_(numNewSamps), adaptive_(adaptive) {
    // Build ROL algorithm and solve SROM optimization problem
    SROMVector<Real> x(Teuchos::rcp(new std::vector<Real>(dimension_*nSamp_,0.)),
                       Teuchos::rcp(new std::vector<Real>(nSamp_,0.)));
    StdVector<Real> l(Teuchos::rcp(new std::vector<Real>(1,0.)));
    bnd_ = Teuchos::rcp(new SROMBoundConstraint<Real>(dimension_));
    con_ = Teuchos::rcp(new SROMEqualityConstraint<Real>);
    bool useAugLag = false;
    buildOptimizer(useAugLag);
    algo_->run(x,l,*obj_,*con_,*bnd_,!SampleGenerator<Real>::batchID());
    // Prune samples with zero weight and set samples/weights
    std::vector<std::vector<Real> > allPoints;
    std::vector<Real> allWeights;
    pruneSamples(allPoints,allWeights,x);
    splitSamples(allPoints,allWeights);
  }

  SROMGenerator(Teuchos::RCP<BatchManager<Real> > &bman, 
                Teuchos::RCP<Objective<Real> > &obj,
                Teuchos::RCP<BoundConstraint<Real> > &bnd,
                const size_t dimension,
                const size_t nSamp = 10,
                const bool adaptive = false,
                const size_t numNewSamps = 0 ) 
    : SampleGenerator<Real>(bman), obj_(obj), bnd_(bnd), dimension_(dimension),
      nSamp_(nSamp), numNewSamps_(numNewSamps), adaptive_(adaptive) {
    // Build ROL algorithm and solve SROM optimization problem
    SROMVector<Real> x(Teuchos::rcp(new std::vector<Real>(dimension_*nSamp_,0.)),
                       Teuchos::rcp(new std::vector<Real>(nSamp_,0.)));
    StdVector<Real> l(Teuchos::rcp(new std::vector<Real>(1,0.)));
    con_ = Teuchos::rcp(new SROMEqualityConstraint<Real>);
    bool useAugLag = false;
    buildOptimizer(useAugLag);
    algo_->run(x,l,*obj_,*con_,*bnd_,!SampleGenerator<Real>::batchID());
    // Prune samples with zero weight and set samples/weights
    std::vector<std::vector<Real> > allPoints;
    std::vector<Real> allWeights;
    pruneSamples(allPoints,allWeights,x);
    splitSamples(allPoints,allWeights);
  }

  SROMGenerator(Teuchos::RCP<BatchManager<Real> > &bman, 
                Teuchos::RCP<Objective<Real> > &obj,
                Teuchos::RCP<BoundConstraint<Real> > &bnd,
                Teuchos::RCP<Vector<Real> > &x,
                const size_t dimension,
                const size_t nSamp = 10,
                const bool adaptive = false,
                const size_t numNewSamps = 0 ) 
    : SampleGenerator<Real>(bman), obj_(obj), bnd_(bnd), dimension_(dimension),
      nSamp_(nSamp), numNewSamps_(numNewSamps), adaptive_(adaptive) {
    // Build ROL algorithm and solve SROM optimization problem
    StdVector<Real> l(Teuchos::rcp(new std::vector<Real>(1,0.)));
    con_ = Teuchos::rcp(new SROMEqualityConstraint<Real>);
    bool useAugLag = true;
    buildOptimizer(useAugLag);
    algo_->run(*x,l,*obj_,*con_,*bnd_,!SampleGenerator<Real>::batchID());
    // Prune samples with zero weight and set samples/weights
    const SROMVector<Real> &ex = Teuchos::dyn_cast<const SROMVector<Real> >(*x);
    std::vector<std::vector<Real> > allPoints;
    std::vector<Real> allWeights;
    pruneSamples(allPoints,allWeights,ex);
    splitSamples(allPoints,allWeights);
  }

  void refine(void) {}
};

}

#endif
