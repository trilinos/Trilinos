// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SROMGENERATOR_HPP
#define ROL_SROMGENERATOR_HPP

#include "ROL_OptimizationSolver.hpp"
#include "ROL_ScalarLinearConstraint.hpp"
#include "ROL_SampleGenerator.hpp"
#include "ROL_MomentObjective.hpp"
#include "ROL_CDFObjective.hpp"
#include "ROL_LinearCombinationObjective.hpp"
#include "ROL_SROMVector.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_SingletonVector.hpp"
#include "ROL_Bounds.hpp"

namespace ROL {

template<class Real>
class SROMGenerator : public SampleGenerator<Real> {
private:
  // Parameterlist for optimization
  ROL::ParameterList parlist_;
  // Vector of distributions (size = dimension of space)
  std::vector<Ptr<Distribution<Real>>> dist_;

  const int dimension_;
  int numSamples_;
  int numMySamples_;
  int numNewSamples_;
  bool adaptive_;

  Real ptol_;
  Real atol_;

  void pruneSamples(const ProbabilityVector<Real> &prob,
                    const AtomVector<Real>        &atom) {
    // Remove points with zero weight
    std::vector<std::vector<Real>> pts;
    std::vector<Real> wts;
    for (int i = 0; i < numMySamples_; i++) {
      if ( prob.getProbability(i) > ptol_ ) {
        pts.push_back(*(atom.getAtom(i)));
        wts.push_back(prob.getProbability(i));
      }
    }
    numMySamples_ = wts.size();
    // Remove atoms that are within atol of each other
    Real err = 0.0;
    std::vector<Real> pt;
    std::vector<int> ind;
    for (int i = 0; i < numMySamples_; i++) {
      pt = pts[i]; ind.clear();
      for (int j = i+1; j < numMySamples_; j++) {
        err = 0.0;
        for (int d = 0; d < dimension_; d++) {
          err += std::pow(pt[d] - pts[j][d],2);
        }
        err = std::sqrt(err);
        if ( err < atol_ ) {
          ind.push_back(j);
          for (int d = 0; d < dimension_; d++) {
            pts[i][d] += pts[j][d];
            wts[i]    += wts[j];
          }
        }
      }
      if ( ind.size() > 0 ) {
        for (int d = 0; d < dimension_; d++) {
          pts[i][d] /= (Real)(ind.size()+1);
        }
        for (int k = ind.size()-1; k >= 0; k--) {
          pts.erase(pts.begin()+ind[k]);
          wts.erase(wts.begin()+ind[k]);
        }
      }
      numMySamples_ = wts.size();
    }
    // Renormalize weights
    Real psum = 0.0, sum = 0.0;
    for (int i = 0; i < numMySamples_; i++) {
      psum += wts[i];
    }
    SampleGenerator<Real>::sumAll(&psum,&sum,1);
    for (int i = 0; i < numMySamples_; i++) {
      wts[i] /= sum;
    }
    // Set points and weights
    SampleGenerator<Real>::setPoints(pts);
    SampleGenerator<Real>::setWeights(wts);
  }

public:

  SROMGenerator(ROL::ParameterList               &parlist,
          const Ptr<BatchManager<Real>>              &bman,
          const std::vector<Ptr<Distribution<Real>>> &dist,
                std::ostream                         &outStream = std::cout)
    : SampleGenerator<Real>(bman), parlist_(parlist), dist_(dist),
      dimension_(dist.size()) {
    // Get SROM sublist
    ROL::ParameterList &list = parlist.sublist("SOL").sublist("Sample Generator").sublist("SROM");
    numSamples_    = list.get("Number of Samples",50);
    adaptive_      = list.get("Adaptive Sampling",false);
    numNewSamples_ = list.get("Number of New Samples Per Adaptation",0);
    ptol_          = list.get("Probability Tolerance",1.e2*std::sqrt(ROL_EPSILON<Real>()));
    atol_          = list.get("Atom Tolerance",1.e2*std::sqrt(ROL_EPSILON<Real>()));
    bool presolve  = list.get("Presolve for Atom Locations",false);
    // Compute batch local number of samples
    int rank  = (int)SampleGenerator<Real>::batchID();
    int nProc = (int)SampleGenerator<Real>::numBatches();
    int frac  = numSamples_ / nProc;
    int rem   = numSamples_ % nProc;
    numMySamples_ = frac + ((rank < rem) ? 1 : 0);
    // Initialize vectors
    Ptr<ProbabilityVector<Real>> prob, prob_lo, prob_hi, prob_eq;
    Ptr<AtomVector<Real>> atom, atom_lo, atom_hi, atom_eq;
    Ptr<Vector<Real>> x, x_lo, x_hi, x_eq;
    initialize_vectors(prob,prob_lo,prob_hi,prob_eq,atom,atom_lo,atom_hi,atom_eq,x,x_lo,x_hi,x_eq,bman);
    Ptr<Vector<Real>> l = makePtr<SingletonVector<Real>>(0.0);
    // Initialize constraints
    Ptr<BoundConstraint<Real>> bnd = makePtr<Bounds<Real>>(x_lo,x_hi);
    Ptr<Constraint<Real>>      con = makePtr<ScalarLinearConstraint<Real>>(x_eq,1.0);
    if (presolve) { // Optimize over atom locations only
      ROL::ParameterList pslist(list);
      pslist.sublist("Step").set("Type","Trust Region");
      Ptr<Objective<Real>> obj = initialize_objective(dist,bman,false,true,pslist);
      OptimizationProblem<Real> optProblem(obj,x,bnd);
      OptimizationSolver<Real>  optSolver(optProblem, pslist);
      optSolver.solve(outStream);
    }
    // Optimization over atom locations and probabilities
    Ptr<Objective<Real>> obj = initialize_objective(dist,bman,true,true,list);
    OptimizationProblem<Real> optProblem(obj,x,bnd,con,l);
    optProblem.check(outStream);
    OptimizationSolver<Real>  optSolver(optProblem, list);
    optSolver.solve(outStream);
    // Prune samples with zero weight and set samples/weights
    pruneSamples(*prob,*atom);
  }

  void refine(void) {}

private:

  void get_scaling_vectors(std::vector<Real> &typw, std::vector<Real> &typx) const {
    typw.clear(); typx.clear();
    typw.resize(numMySamples_,(Real)(numSamples_*numSamples_));
    typx.resize(numMySamples_*dimension_,0);
    Real mean = 1, var = 1, one(1);
    for (int j = 0; j < dimension_; j++) {
      mean = std::abs(dist_[j]->moment(1));
      var  = dist_[j]->moment(2) - mean*mean;
      mean = ((mean > ROL_EPSILON<Real>()) ? mean : std::sqrt(var));
      mean = ((mean > ROL_EPSILON<Real>()) ? mean : one);
      for (int i = 0; i < numMySamples_; i++) {
        typx[i*dimension_ + j] = one/(mean*mean);
      }
    }
  }

  void initialize_vectors(Ptr<ProbabilityVector<Real>>  &prob,
                          Ptr<ProbabilityVector<Real>>  &prob_lo,
                          Ptr<ProbabilityVector<Real>>  &prob_hi,
                          Ptr<ProbabilityVector<Real>>  &prob_eq,
                          Ptr<AtomVector<Real>>         &atom,
                          Ptr<AtomVector<Real>>         &atom_lo,
                          Ptr<AtomVector<Real>>         &atom_hi,
                          Ptr<AtomVector<Real>>         &atom_eq,
                          Ptr<Vector<Real>>             &vec,
                          Ptr<Vector<Real>>             &vec_lo,
                          Ptr<Vector<Real>>             &vec_hi,
                          Ptr<Vector<Real>>             &vec_eq,
                          const Ptr<BatchManager<Real>> &bman) const {
    // Compute scaling for probability and atom vectors
    std::vector<Real> typx, typw;
    get_scaling_vectors(typw,typx);
    // Compute initial guess and bounds for probability and atom vectors
    std::vector<Real> pt(dimension_*numMySamples_,0.), wt(numMySamples_,1./(Real)numSamples_);
    std::vector<Real> pt_lo(dimension_*numMySamples_,0.), pt_hi(dimension_*numMySamples_,0.);
    std::vector<Real> wt_lo(numMySamples_,0.), wt_hi(numMySamples_,1.);
    std::vector<Real> pt_eq(dimension_*numMySamples_,0.), wt_eq(numMySamples_,1.);
    Real lo = 0., hi = 0.;
    srand(12345*SampleGenerator<Real>::batchID());
    for ( int j = 0; j < dimension_; j++) {
      lo = dist_[j]->lowerBound();
      hi = dist_[j]->upperBound();
      for (int i = 0; i < numMySamples_; i++) {
        pt[i*dimension_ + j] = dist_[j]->invertCDF((Real)rand()/(Real)RAND_MAX);
        //pt[i*dimension_ + j] = dist_[j]->invertCDF(0);
        pt_lo[i*dimension_ + j] = lo;
        pt_hi[i*dimension_ + j] = hi;
      }
    }
    // Build probability, atom, and SROM vectors 
    prob = makePtr<PrimalProbabilityVector<Real>>(
           makePtr<std::vector<Real>>(wt),bman,
           makePtr<std::vector<Real>>(typw));
    atom = makePtr<PrimalAtomVector<Real>>(
           makePtr<std::vector<Real>>(pt),bman,numMySamples_,dimension_,
           makePtr<std::vector<Real>>(typx));
    vec  = makePtr<SROMVector<Real>>(prob,atom);
    // Lower and upper bounds on Probability Vector
    prob_lo = makePtr<PrimalProbabilityVector<Real>>(
              makePtr<std::vector<Real>>(wt_lo),bman,
              makePtr<std::vector<Real>>(typw));
    prob_hi = makePtr<PrimalProbabilityVector<Real>>(
              makePtr<std::vector<Real>>(wt_hi),bman,
              makePtr<std::vector<Real>>(typw));
    // Lower and upper bounds on Atom Vector
    atom_lo = makePtr<PrimalAtomVector<Real>>(
              makePtr<std::vector<Real>>(pt_lo),bman,numMySamples_,dimension_,
              makePtr<std::vector<Real>>(typx));
    atom_hi = makePtr<PrimalAtomVector<Real>>(
              makePtr<std::vector<Real>>(pt_hi),bman,numMySamples_,dimension_,
              makePtr<std::vector<Real>>(typx));
   // Lower and upper bounds on SROM Vector
    vec_lo = makePtr<SROMVector<Real>>(prob_lo,atom_lo);
    vec_hi = makePtr<SROMVector<Real>>(prob_hi,atom_hi);
    // Constraint vectors
    prob_eq = makePtr<DualProbabilityVector<Real>>(
              makePtr<std::vector<Real>>(wt_eq),bman,
              makePtr<std::vector<Real>>(typw));
    atom_eq = makePtr<DualAtomVector<Real>>(
              makePtr<std::vector<Real>>(pt_eq),bman,numMySamples_,dimension_,
              makePtr<std::vector<Real>>(typx));
    vec_eq  = makePtr<SROMVector<Real>>(prob_eq,atom_eq);
  }

  Ptr<Objective<Real>>  initialize_objective(const std::vector<Ptr<Distribution<Real>>> &dist,
                                             const Ptr<BatchManager<Real>>              &bman,
                                             const bool optProb, const bool optAtom,
                                             ROL::ParameterList                     &list) const {
    std::vector<Ptr<Objective<Real>>> obj_vec;
    // Build CDF objective function
    Real scale = list.get("CDF Smoothing Parameter",1.e-2);
    obj_vec.push_back(makePtr<CDFObjective<Real>>(dist,bman,scale,optProb,optAtom));
    // Build moment matching objective function
    std::vector<int> tmp_order
      = ROL::getArrayFromStringParameter<int>(list,"Moments");
    std::vector<int> order(tmp_order.size(),0);
    for (unsigned int i = 0; i < tmp_order.size(); i++) {
      order[i] = static_cast<int>(tmp_order[i]);
    }
    obj_vec.push_back(makePtr<MomentObjective<Real>>(dist,order,bman,optProb,optAtom));
    // Build linear combination objective function
    std::vector<Real> tmp_coeff
      = ROL::getArrayFromStringParameter<Real>(list,"Coefficients");
    std::vector<Real> coeff(2,0.);
    coeff[0] = tmp_coeff[0]; coeff[1] = tmp_coeff[1];
    return makePtr<LinearCombinationObjective<Real>>(coeff,obj_vec);
  }

  int numGlobalSamples(void) const {
    return numSamples_;
  }
};

}

#endif
