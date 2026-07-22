// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_OED_GREEDYALGORITHM_HPP
#define ROL_OED_GREEDYALGORITHM_HPP

#include <iostream>

#include "ROL_Ptr.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_LinearAlgebra.hpp"
#include "ROL_BLAS.hpp"
#include "ROL_LAPACK.hpp"

#include "ROL_Vector.hpp"
#include "ROL_Objective.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_SampleGenerator.hpp"
#include "ROL_VectorController.hpp"
#include "ROL_ObjectiveFromConstraint.hpp"

#include "ROL_OED_Factors.hpp"
#include "ROL_OED_MomentOperator.hpp"
#include "ROL_OED_TraceSampler.hpp"
#include "ROL_OED_Radamacher.hpp"

namespace ROL::OED {

template<typename Real>
struct GreedyObjective {
public:
  const Ptr<MomentOperator<Real>> cov;
  const Ptr<Factors<Real>> factors;
  const Ptr<BatchManager<Real>> bman;
  const bool useDeletion;

  GreedyObjective(const Ptr<MomentOperator<Real>>& cov_, bool useDeletion_)
    : cov(cov_), factors(cov_->getFactors()), useDeletion(useDeletion_),
      lapack(makePtr<LAPACK<int,Real>>()), blas(makePtr<BLAS<int,Real>>()),
      nobs(cov->getFactors()->numObservations()) {
    obs.resize(nobs);
    for (unsigned i=0u; i<nobs; ++i) {
      obs[i] = factors->createObservationVector(true);
      obs[i]->set(*obs[i]->basis(i));
    }
  }
  virtual ~GreedyObjective() {}
  virtual void update(const Vector<Real>& x, int iter=-1) {
    cov->update(x,UpdateType::Trial,iter);
  }
  virtual Real value(const Vector<Real>& x) const = 0;
  virtual Real computeIncrement(const Vector<Real>& x, unsigned ind) const = 0;
  virtual std::string getName() const { return "Undefined"; }

protected:
  const Ptr<LAPACK<int,Real>> lapack;
  const Ptr<BLAS<int,Real>> blas;
  const unsigned nobs;
  std::vector<Ptr<Vector<Real>>> obs;
};

template<typename Real>
struct GreedyObjectiveA : public GreedyObjective<Real> {
public:
  GreedyObjectiveA(const Ptr<MomentOperator<Real>>& cov_, const Ptr<TraceSampler<Real>>& tsampler, std::vector<Real> weight, bool useDeletion_)
    : GreedyObjective<Real>(cov_,useDeletion_), tsampler_(tsampler),
      weight_(weight), size_(weight.size()),
      stateStore_(makePtr<VectorController<Real,unsigned>>()),
      uvec_(cov_->getFactors()->createParameterVector(false)),
      rhs_(uvec_->dual().clone()) {
    X_.resize(nobs);
    Y_.resize(nobs);
    XY_.reshape(nobs,nobs);
    YY_.reshape(nobs,nobs);
    for (unsigned i=0; i<nobs; ++i) {
      X_[i] = GreedyObjective<Real>::factors->createParameterVector(true);
      Y_[i] = GreedyObjective<Real>::factors->createParameterVector(false);
    }
  }
  void update(const Vector<Real>& x, int iter=-1) override final {
    GreedyObjective<Real>::update(x,iter);
    stateStore_->objectiveUpdate(UpdateType::Trial);
    val_ = static_cast<Real>(0);
    for (unsigned i = 0; i < size_; ++i) {
      tsampler_->get(*rhs_,{static_cast<Real>(i)});
      GreedyObjective<Real>::cov->applyInverse(*uvec_,*rhs_,x);
      stateStore_->set(*uvec_,i);
      val_ += weight_[i] * rhs_->apply(*uvec_);
    }
  }
  Real value(const Vector<Real>& x) const override final {
    return val_;
  }
  Real computeIncrement(const Vector<Real>& x, unsigned ind) const override final {
    const Real one(1);
    const Real scal(GreedyObjective<Real>::useDeletion ? -one : one);
    for (unsigned i=0; i<nobs; ++i) {
      GreedyObjective<Real>::factors->applyAdjoint(*X_[i],*obs[i],ind); // Can Compute X offline
      GreedyObjective<Real>::cov->applyInverse(*Y_[i],*X_[i],x);
      XY_(i,i) = one + scal*X_[i]->apply(*Y_[i]);
      YY_(i,i) = Y_[i]->dot(*Y_[i]);
      for (unsigned j=0u; j<i; ++j) {
        XY_(i,j) = scal*X_[i]->apply(*Y_[j]);
        YY_(i,j) = Y_[i]->dot(*Y_[j]);
        YY_(j,i) = YY_(i,j);
      }
    }
    int info;
    lapack->POTRF('L',nobs,XY_.values(),nobs,&info);
    lapack->POTRS('L',nobs,nobs,XY_.values(),nobs,YY_.values(),nobs,&info);
    Real inc(0);
    for (unsigned i=0u; i<nobs; ++i)
      inc += YY_(i,i);
    return inc;
  }
  std::string getName() const override final { return "A-Optimality"; }
private:
  const Ptr<TraceSampler<Real>> tsampler_;
  const std::vector<Real> weight_;
  const unsigned size_;
  const Ptr<VectorController<Real,unsigned>> stateStore_;
  const Ptr<Vector<Real>> uvec_, rhs_;
  Real val_;
  std::vector<Ptr<Vector<Real>>> X_, Y_;
  mutable LA::Matrix<Real> XY_, YY_;

  using GreedyObjective<Real>::lapack;
  using GreedyObjective<Real>::nobs;
  using GreedyObjective<Real>::obs;
};

template<typename Real>
struct GreedyObjectiveC : public GreedyObjective<Real> {
public:
  GreedyObjectiveC(const Ptr<MomentOperator<Real>>& cov_, const Ptr<Vector<Real>>& c, bool useDeletion_)
    : GreedyObjective<Real>(cov_,useDeletion_), c_(c),
      uvec_(cov_->getFactors()->createParameterVector(false)),
      Y_(cov_->getFactors()->createParameterVector(false)) {
    X_.resize(nobs);
    XY_.reshape(nobs,nobs);
    Yc_.resize(nobs);
    SYc_.resize(nobs);
    for (unsigned i=0; i<nobs; ++i)
      X_[i] = GreedyObjective<Real>::factors->createParameterVector(true);
  }
  void update(const Vector<Real>& x, int iter=-1) {
    GreedyObjective<Real>::update(x,iter);
    GreedyObjective<Real>::cov->applyInverse(*uvec_,*c_,x);
  }
  Real value(const Vector<Real>& x) const override final {
    return c_->apply(*uvec_);
  }
  Real computeIncrement(const Vector<Real>& x, unsigned ind) const override final {
    const Real one(1);
    const Real scal(GreedyObjective<Real>::useDeletion ? -one : one);
    for (unsigned i=0; i<nobs; ++i) {
      GreedyObjective<Real>::factors->applyAdjoint(*X_[i],*obs[i],ind); // Can Compute X offline
      GreedyObjective<Real>::cov->applyInverse(*Y_,*X_[i],x);
      Yc_(i) = c_->apply(*Y_);
      SYc_(i) = Yc_(i);
      for (unsigned j=0u; j<=i; ++j)
        XY_(i,j) = scal*Y_->apply(*X_[j]);
      XY_(i,i) += one;
    }
    int info;
    lapack->POTRF('L',nobs,XY_.values(),nobs,&info);
    lapack->POTRS('L',nobs,1,XY_.values(),nobs,SYc_.values(),nobs,&info);
    return blas->DOT(nobs,Yc_.values(),1,SYc_.values(),1);
  }
  std::string getName() const override final { return "C-Optimality"; }
private:
  const Ptr<Vector<Real>> c_, uvec_, Y_;
  std::vector<Ptr<Vector<Real>>> X_;
  mutable LA::Matrix<Real> XY_;
  mutable LA::Vector<Real> Yc_, SYc_;

  using GreedyObjective<Real>::lapack;
  using GreedyObjective<Real>::blas;
  using GreedyObjective<Real>::nobs;
  using GreedyObjective<Real>::obs;
};

template<typename Real>
struct GreedyObjectiveD : public GreedyObjective<Real> {
public:
  GreedyObjectiveD(const Ptr<MomentOperator<Real>>& cov_, bool useDeletion_)
    : GreedyObjective<Real>(cov_,useDeletion_),
      Y_(cov_->getFactors()->createParameterVector(false)) {
    X_.resize(nobs);
    XY_.reshape(nobs,nobs);
    for (unsigned i=0; i<nobs; ++i)
      X_[i] = GreedyObjective<Real>::factors->createParameterVector(true);
  }
  Real value(const Vector<Real>& x) const override final {
    return -GreedyObjective<Real>::cov->logDeterminant(x);
  }
  Real computeIncrement(const Vector<Real>& x, unsigned ind) const override final {
    const Real one(1);
    const Real scal(GreedyObjective<Real>::useDeletion ? -one : one);
    for (unsigned i=0; i<nobs; ++i) {
      GreedyObjective<Real>::factors->applyAdjoint(*X_[i],*obs[i],ind); // Can Compute X offline
      GreedyObjective<Real>::cov->applyInverse(*Y_,*X_[i],x);
      for (unsigned j=0u; j<=i; ++j)
        XY_(i,j) = scal*Y_->apply(*X_[j]);
      XY_(i,i) += one;
    }
    int info;
    lapack->POTRF('L',nobs,XY_.values(),nobs,&info);
    Real inc(1);
    for (unsigned i=0u; i<nobs; ++i)
      inc *= XY_(i,i)*XY_(i,i);
    return -scal*(one-inc);
  }
  std::string getName() const override final { return "D-Optimality"; }
private:
  const Ptr<Vector<Real>> Y_;
  std::vector<Ptr<Vector<Real>>> X_;
  mutable LA::Matrix<Real> XY_;

  using GreedyObjective<Real>::lapack;
  using GreedyObjective<Real>::nobs;
  using GreedyObjective<Real>::obs;
};

template<typename Real>
struct GreedyObjectiveI : public GreedyObjective<Real> {
public:
  GreedyObjectiveI(const Ptr<MomentOperator<Real>>& cov_,
                   const Ptr<Factors<Real>>& PVfactors, const Ptr<SampleGenerator<Real>>& sampler,
                   const Ptr<TraceSampler<Real>>& tsampler, std::vector<Real> weight,
                         bool useDeletion_)
    : GreedyObjective<Real>(cov_,useDeletion_), tsampler_(tsampler),
      size_(weight.size()),
      uStateStore_(makePtr<VectorController<Real,unsigned>>()),
      vStateStore_(makePtr<VectorController<Real,unsigned>>()),
      uvec_(cov_->getFactors()->getTheta()->clone()),
      vvec_(uvec_->clone()), rhs_(uvec_->dual().clone()),
      Y_(uvec_->clone()) {
    auto F = PVfactors->createParameterVector(true);
    auto b = F->clone();
    auto obs0 = PVfactors->createObservationVector(false);
    b_.clear(); b_.resize(size_);
    for (unsigned i=0u; i<size_; ++i) {
      b->zero();
      b_[i] = b->clone();
      tsampler->get(*rhs_,{static_cast<Real>(i)});
      for (int j = 0; j < sampler->numMySamples(); ++j) {
        PVfactors->apply(*obs0,rhs_->dual(),sampler->getMyPoint(j));
	PVfactors->applyAdjoint(*F,obs0->dual(),sampler->getMyPoint(j));
        b->axpy(weight[i]*sampler->getMyWeight(j),*F);
      }
      sampler->sumAll(*b,*b_[i]);
    }
    X_.resize(nobs);
    XY_.reshape(nobs,nobs);
    Xu_.resize(nobs);
    Xv_.resize(nobs);
    for (unsigned i=0; i<nobs; ++i)
      X_[i] = GreedyObjective<Real>::factors->createParameterVector(true);
  }
  void update(const Vector<Real>& x, int iter=-1) override final {
    GreedyObjective<Real>::update(x,iter);
    uStateStore_->objectiveUpdate(UpdateType::Trial);
    vStateStore_->objectiveUpdate(UpdateType::Trial);
    val_ = static_cast<Real>(0);
    for (unsigned i = 0; i < size_; ++i) {
      tsampler_->get(*rhs_,{static_cast<Real>(i)});
      GreedyObjective<Real>::cov->applyInverse(*uvec_,*rhs_,x);
      uStateStore_->set(*uvec_,i);
      val_ += b_[i]->apply(*uvec_);
      GreedyObjective<Real>::cov->applyInverse(*vvec_,*b_[i],x);
      vStateStore_->set(*vvec_,i);
    }
  }
  Real value(const Vector<Real>& x) const override final {
    return val_;
  }
  Real computeIncrement(const Vector<Real>& x, unsigned ind) const override final {
    const Real one(1);
    const Real scal(GreedyObjective<Real>::useDeletion ? -one : one);
    for (unsigned i=0u; i<nobs; ++i) {
      GreedyObjective<Real>::factors->applyAdjoint(*X_[i],*obs[i],ind); // Can Compute X offline
      GreedyObjective<Real>::cov->applyInverse(*Y_,*X_[i],x);
      for (unsigned j=0u; j<=i; ++j)
        XY_(i,j) = scal*Y_->apply(*X_[j]);
      XY_(i,i) += one;
    }
    int infoF, infoS;
    lapack->POTRF('L',nobs,XY_.values(),nobs,&infoF);
    Real inc(0);
    for (unsigned i=0u; i<size_; ++i) {
      uStateStore_->get(*uvec_,i);
      vStateStore_->get(*vvec_,i);
      for (unsigned j=0u; j<nobs; ++j) {
        Xu_(j) = X_[j]->apply(*uvec_);
        Xv_(j) = X_[j]->apply(*vvec_);
      }
      lapack->POTRS('L',nobs,1,XY_.values(),nobs,Xu_.values(),nobs,&infoS);
      inc += blas->DOT(nobs,Xu_.values(),1,Xv_.values(),1);
    }
    return inc;
  }
  std::string getName() const override final { return "I-Optimality"; }
private:
  const Ptr<TraceSampler<Real>> tsampler_;
  const unsigned size_;
  const Ptr<VectorController<Real,unsigned>> uStateStore_, vStateStore_;
  const Ptr<Vector<Real>> uvec_, vvec_, rhs_, Y_;
  std::vector<Ptr<Vector<Real>>> b_, X_;
  Real val_;
  mutable LA::Matrix<Real> XY_;
  mutable LA::Vector<Real> Xv_, Xu_;

  using GreedyObjective<Real>::lapack;
  using GreedyObjective<Real>::blas;
  using GreedyObjective<Real>::nobs;
  using GreedyObjective<Real>::obs;
};

// General Greedy Algorithm Implementation
template<typename Real>
Ptr<Vector<Real>> GreedyAlgorithm(const Ptr<GreedyObjective<Real>>& obj,
                                  const Ptr<const Vector<Real>>& cvec,
                                        Real budget,
                                        ParameterList& list,
                                        std::ostream& os=std::cout) {
  // Set up algorithm for either adding or deleting experiments
  bool useDeletion = list.sublist("OED").sublist("Greedy").get("Use Deletion",true);
  Real scale = useDeletion ? static_cast<Real>(-1) : static_cast<Real>(1);
  Real vval  = useDeletion ? static_cast<Real>(0)  : static_cast<Real>(1);

  // Set up solution and cost vectors
  auto pvec = cvec->dual().clone();
  pvec->setScalar(static_cast<Real>(1)-vval);
  const auto& pdata = staticPtrCast<ProbabilityVector<Real>>(pvec)->getVector();
  const auto& cdata = staticPtrCast<const ProbabilityVector<Real>>(cvec)->getVector();

  // Get parallel processing information
  const auto& bman = staticPtrCast<ProbabilityVector<Real>>(pvec)->getBatchManager();
  const unsigned nbatch = bman->numBatches();
  const unsigned rank   = bman->batchID();
  const unsigned size   = cvec->dimension();

  // Set up index vector
  const unsigned mysize = staticPtrCast<ProbabilityVector<Real>>(pvec)->getNumMyAtoms();
  std::vector<unsigned> myind(mysize);
  for (unsigned i=0u; i<mysize; ++i) myind[i] = i;

  // Print output to stream
  auto type = obj->getName();
  os << std::endl;
  if (useDeletion)
    os << "Greedy Deletion Algorithm for Binary " << type << std::endl;
  else
    os << "Greedy Addition Algorithm for Binary " << type << std::endl;
  os << "  ";
  os << std::setw(6)  << std::left << "iter";
  os << std::setw(15) << std::left << "value";
  os << std::setw(15) << std::left << "cost";
  os << std::setw(15) << std::left << "delta";
  os << std::setw(10) << std::left << "index";
  if (nbatch > 1u)
    os << std::setw(10) << std::left << "rank";
  os << std::endl;

  // Run greedy deletion algorithm
  const Real zero(0);
  unsigned jind(0), i(0), mind(0), rind(0);
  Real val(0), vi(0), vij(0);
  Real cost(cvec->apply(*pvec)), cdel(0), ctmp(0);
  std::vector<Real> gvi(4*nbatch), tmp(4);
  for (i=0u; i<size; ++i) {
    obj->update(*pvec,i);
    val = obj->value(*pvec);

    // Print output to stream
    os << "  ";
    os << std::setw(6)  << std::left << i;
    os << std::setw(15) << std::left << val;
    os << std::setw(15) << std::left << cost;
    if (i==0u) {
      os << std::setw(15) << std::left << "---";
      os << std::setw(10) << std::left << "---";
      if (nbatch > 1u)
        os << std::setw(10) << std::left << "---";
    }
    else {
      os << std::setw(15) << std::left << vi;
      os << std::setw(10) << std::left << mind;
      if (nbatch > 1u)
        os << std::setw(10) << std::left << rind;
    }
    os << std::endl;

    // Compute effectivity index
    vi = scale*ROL_NINF<Real>();
    jind = mysize;
    for (unsigned j=0u; j<myind.size(); ++j) {
      vij = obj->computeIncrement(*pvec,myind[j]) / (*cdata)[myind[j]];
      if (j==0u || scale*vi < scale*vij) {
        vi = vij;
        jind = j;
      }
    }

    // Communicate index accross processors
    bman->barrier();
    tmp[0] = vi;
    tmp[1] = static_cast<Real>(jind);
    tmp[2] = static_cast<Real>(myind[jind]);
    tmp[3] = static_cast<Real>(rank);
    bman->gather(tmp.data(),4,gvi.data(),4*nbatch,0);
    if (rank==0u) {
      tmp[0] = gvi[0];
      tmp[1] = gvi[1];
      tmp[2] = gvi[2];
      tmp[3] = gvi[3];
      for (unsigned j=1u; j<nbatch; ++j) {
        if (scale*tmp[0] < scale*gvi[4*j]) {
          tmp[0] = gvi[4*j];
          tmp[1] = gvi[4*j+1];
          tmp[2] = gvi[4*j+2];
          tmp[3] = gvi[4*j+3];
        }
      }
    }
    bman->broadcast(tmp.data(),4,0);
    vi = tmp[0];
    jind = static_cast<unsigned>(tmp[1]);
    mind = static_cast<unsigned>(tmp[2]);
    rind = static_cast<unsigned>(tmp[3]);

    // Update design and cost
    ctmp = zero; cdel = zero;
    if (rank==rind) ctmp = (*cdata)[mind];
    bman->sumAll(&ctmp,&cdel,1);
    if (!useDeletion && cost+cdel > budget) break;
    if (rank==rind) {
      (*pdata)[mind] = static_cast<Real>(vval); // Remove experiment
      myind.erase(myind.begin()+jind);          // Remove experiment
    }
    cost += scale*cdel;                         // Update cost
    if (useDeletion && cost <= budget) break;
  }

  // Print output to stream
  if (useDeletion) {
    i++;
    obj->update(*pvec,i);
    val = obj->value(*pvec);
    os << "  ";
    os << std::setw(6)  << std::left << i;
    os << std::setw(15) << std::left << val;
    os << std::setw(15) << std::left << cost;
    os << std::setw(15) << std::left << vi;
    os << std::setw(10) << std::left << mind;
    if (nbatch > 1u)
      os << std::setw(10) << std::left << rind;
    os << std::endl;
  }
  return pvec;
}

// Specialized Greedy Algorithm for A-Optimality
template<typename Real>
Ptr<Vector<Real>> GreedyAlgorithmA(const Ptr<Constraint<Real>>& model,
                                   const Ptr<Vector<Real>>& obs,
                                   const Ptr<SampleGenerator<Real>>& sampler,
                                   const Ptr<Vector<Real>>& theta,
                                   const Ptr<MomentOperator<Real>>& cov,
                                   const Ptr<Vector<Real>>& cvec,
                                         Real budget,
                                         ParameterList& list,
                                         std::ostream& os=std::cout) {
  list.sublist("OED").set("Optimality Type","A");
  cov->generateFactors(model,theta,obs,sampler);
  bool useRandomTrace = list.sublist("OED").sublist("A-Optimality").get("Randomized Trace Estimation",false);
  int nRandomTrace    = list.sublist("OED").sublist("A-Optimality").get("Number of Samples",100);
  int nfactors        = theta->dimension();
  int size            = (useRandomTrace ? nRandomTrace : nfactors);
  Ptr<TraceSampler<Real>> traceSampler;
  if (useRandomTrace) traceSampler = makePtr<Radamacher<Real>>(theta,size);
  else                traceSampler = makePtr<TraceSampler<Real>>(theta);
  const int one(1);
  Real val = (useRandomTrace ? one/static_cast<Real>(nRandomTrace) : one);
  std::vector<Real> weight(size,val);
  bool useDeletion = list.sublist("OED").sublist("Greedy").get("Use Deletion",true);
  auto obj = makePtr<GreedyObjectiveA<Real>>(cov,traceSampler,weight,useDeletion);
  return GreedyAlgorithm<Real>(obj,cvec,budget,list,os);
}

// Specialized Greedy Algorithm for A-Optimality
template<typename Real>
Ptr<Vector<Real>> GreedyAlgorithmA(const Ptr<Objective<Real>>& model,
                                   const Ptr<SampleGenerator<Real>>& sampler,
                                   const Ptr<Vector<Real>>& theta,
                                   const Ptr<MomentOperator<Real>>& cov,
                                   const Ptr<Vector<Real>>& cvec,
                                         Real budget,
                                         ParameterList& list,
                                         std::ostream& os=std::cout) {
  list.sublist("OED").set("Optimality Type","A");
  cov->generateFactors(model,theta,sampler);
  bool useRandomTrace = list.sublist("OED").sublist("A-Optimality").get("Randomized Trace Estimation",false);
  int nRandomTrace    = list.sublist("OED").sublist("A-Optimality").get("Number of Samples",100);
  int nfactors        = theta->dimension();
  int size            = (useRandomTrace ? nRandomTrace : nfactors);
  Ptr<TraceSampler<Real>> traceSampler;
  if (useRandomTrace) traceSampler = makePtr<Radamacher<Real>>(theta,size);
  else                traceSampler = makePtr<TraceSampler<Real>>(theta);
  const int one(1);
  Real val = (useRandomTrace ? one/static_cast<Real>(nRandomTrace) : one);
  std::vector<Real> weight(size,val);
  bool useDeletion = list.sublist("OED").sublist("Greedy").get("Use Deletion",true);
  auto obj = makePtr<GreedyObjectiveA<Real>>(cov,traceSampler,weight,useDeletion);
  return GreedyAlgorithm<Real>(obj,cvec,budget,list,os);
}

// Specialized Greedy Algorithm for C-Optimality
template<typename Real>
Ptr<Vector<Real>> GreedyAlgorithmC(const Ptr<Constraint<Real>>& model,
                                   const Ptr<Vector<Real>>& obs,
                                   const Ptr<SampleGenerator<Real>>& sampler,
                                   const Ptr<Vector<Real>>& theta,
                                   const Ptr<MomentOperator<Real>>& cov,
                                   const Ptr<Vector<Real>>& CoptVec,
                                   const Ptr<Vector<Real>>& cvec,
                                         Real budget,
                                         ParameterList& list,
                                         std::ostream& os=std::cout) {
  list.sublist("OED").set("Optimality Type","C");
  cov->generateFactors(model,theta,obs,sampler);
  bool useDeletion = list.sublist("OED").sublist("Greedy").get("Use Deletion",true);
  auto obj = makePtr<GreedyObjectiveC<Real>>(cov,CoptVec,useDeletion);
  return GreedyAlgorithm<Real>(obj,cvec,budget,list,os);
}

// Specialized Greedy Algorithm for C-Optimality
template<typename Real>
Ptr<Vector<Real>> GreedyAlgorithmC(const Ptr<Constraint<Real>>& model,
                                   const Ptr<Vector<Real>>& obs,
                                   const Ptr<SampleGenerator<Real>>& sampler,
                                   const Ptr<Vector<Real>>& theta,
                                   const Ptr<MomentOperator<Real>>& cov,
                                   const Ptr<Vector<Real>>& cvec,
                                         Real budget,
                                         ParameterList& list,
                                         std::ostream& os=std::cout) {
  auto c = theta->dual().clone();
  c->setScalar(static_cast<Real>(1));
  return GreedyAlgorithmC<Real>(model,obs,theta,cov,c,cvec,budget,list,os);
}

// Specialized Greedy Algorithm for C-Optimality
template<typename Real>
Ptr<Vector<Real>> GreedyAlgorithmC(const Ptr<Objective<Real>>& model,
                                   const Ptr<SampleGenerator<Real>>& sampler,
                                   const Ptr<Vector<Real>>& theta,
                                   const Ptr<MomentOperator<Real>>& cov,
                                   const Ptr<Vector<Real>>& CoptVec,
                                   const Ptr<Vector<Real>>& cvec,
                                         Real budget,
                                         ParameterList& list,
                                         std::ostream& os=std::cout) {
  list.sublist("OED").set("Optimality Type","C");
  cov->generateFactors(model,theta,sampler);
  bool useDeletion = list.sublist("OED").sublist("Greedy").get("Use Deletion",true);
  auto obj = makePtr<GreedyObjectiveC<Real>>(cov,CoptVec,useDeletion);
  return GreedyAlgorithm<Real>(obj,cvec,budget,list,os);
}

// Specialized Greedy Algorithm for C-Optimality
template<typename Real>
Ptr<Vector<Real>> GreedyAlgorithmC(const Ptr<Objective<Real>>& model,
                                   const Ptr<SampleGenerator<Real>>& sampler,
                                   const Ptr<Vector<Real>>& theta,
                                   const Ptr<MomentOperator<Real>>& cov,
                                   const Ptr<Vector<Real>>& cvec,
                                         Real budget,
                                         ParameterList& list,
                                         std::ostream& os=std::cout) {
  auto c = theta->dual().clone();
  c->setScalar(static_cast<Real>(1));
  return GreedyAlgorithmC<Real>(model,theta,cov,c,cvec,budget,list,os);
}

// Specialized Greedy Algorithm for D-Optimality
template<typename Real>
Ptr<Vector<Real>> GreedyAlgorithmD(const Ptr<Constraint<Real>>& model,
                                   const Ptr<Vector<Real>>& obs,
                                   const Ptr<SampleGenerator<Real>>& sampler,
                                   const Ptr<Vector<Real>>& theta,
                                   const Ptr<MomentOperator<Real>>& cov,
                                   const Ptr<Vector<Real>>& cvec,
                                         Real budget,
                                         ParameterList& list,
                                         std::ostream& os=std::cout) {
  list.sublist("OED").set("Optimality Type","D");
  cov->generateFactors(model,theta,obs,sampler);
  bool useDeletion = list.sublist("OED").sublist("Greedy").get("Use Deletion",true);
  auto obj = makePtr<GreedyObjectiveD<Real>>(cov,useDeletion);
  return GreedyAlgorithm<Real>(obj,cvec,budget,list,os);
}

// Specialized Greedy Algorithm for D-Optimality
template<typename Real>
Ptr<Vector<Real>> GreedyAlgorithmD(const Ptr<Objective<Real>>& model,
                                   const Ptr<SampleGenerator<Real>>& sampler,
                                   const Ptr<Vector<Real>>& theta,
                                   const Ptr<MomentOperator<Real>>& cov,
                                   const Ptr<Vector<Real>>& cvec,
                                         Real budget,
                                         ParameterList& list,
                                         std::ostream& os=std::cout) {
  list.sublist("OED").set("Optimality Type","D");
  cov->generateFactors(model,theta,sampler);
  bool useDeletion = list.sublist("OED").sublist("Greedy").get("Use Deletion",true);
  auto obj = makePtr<GreedyObjectiveD<Real>>(cov,useDeletion);
  return GreedyAlgorithm<Real>(obj,cvec,budget,list,os);
}

// Specialized Greedy Algorithm for I-Optimality
template<typename Real>
Ptr<Vector<Real>> GreedyAlgorithmI(const Ptr<Constraint<Real>>& model,
                                   const Ptr<Vector<Real>>& obs,
                                   const Ptr<SampleGenerator<Real>>& sampler,
                                   const Ptr<Vector<Real>>& theta,
                                   const Ptr<MomentOperator<Real>>& cov,
                                   const Ptr<Constraint<Real>>& Ifun,
                                   const Ptr<Vector<Real>>& Ivec,
                                   const Ptr<SampleGenerator<Real>>& Isampler,
                                   const Ptr<Vector<Real>>& cvec,
                                         Real budget,
                                         ParameterList& list,
                                         std::ostream& os=std::cout) {
  list.sublist("OED").set("Optimality Type","I");
  cov->generateFactors(model,theta,obs,sampler);
  auto factorsPV = makePtr<Factors<Real>>(Ifun,theta,Ivec,Isampler);
  bool useRandomTrace = list.sublist("OED").sublist("I-Optimality").get("Randomized Trace Estimation",false);
  int nRandomTrace    = list.sublist("OED").sublist("I-Optimality").get("Number of Samples",100);
  int nfactors        = theta->dimension();
  int size            = (useRandomTrace ? nRandomTrace : nfactors);
  Ptr<TraceSampler<Real>> traceSampler;
  if (useRandomTrace) traceSampler = makePtr<Radamacher<Real>>(theta,size);
  else                traceSampler = makePtr<TraceSampler<Real>>(theta);
  const int one(1);
  Real val = (useRandomTrace ? one/static_cast<Real>(nRandomTrace) : one);
  std::vector<Real> weight(size,val);
  bool useDeletion = list.sublist("OED").sublist("Greedy").get("Use Deletion",true);
  auto obj = makePtr<GreedyObjectiveI<Real>>(cov,factorsPV,Isampler,traceSampler,weight,useDeletion);
  return GreedyAlgorithm<Real>(obj,cvec,budget,list,os);
}

// Specialized Greedy Algorithm for I-Optimality
template<typename Real>
Ptr<Vector<Real>> GreedyAlgorithmI(const Ptr<Objective<Real>>& model,
                                   const Ptr<SampleGenerator<Real>>& sampler,
                                   const Ptr<Vector<Real>>& theta,
                                   const Ptr<MomentOperator<Real>>& cov,
                                   const Ptr<Constraint<Real>>& Ifun,
                                   const Ptr<Vector<Real>>& Ivec,
                                   const Ptr<SampleGenerator<Real>>& Isampler,
                                   const Ptr<Vector<Real>>& cvec,
                                         Real budget,
                                         ParameterList& list,
                                         std::ostream& os=std::cout) {
  list.sublist("OED").set("Optimality Type","I");
  cov->generateFactors(model,theta,sampler);
  auto factorsPV = makePtr<Factors<Real>>(Ifun,theta,Ivec,Isampler);
  bool useRandomTrace = list.sublist("OED").sublist("I-Optimality").get("Randomized Trace Estimation",false);
  int nRandomTrace    = list.sublist("OED").sublist("I-Optimality").get("Number of Samples",100);
  int nfactors        = theta->dimension();
  int size            = (useRandomTrace ? nRandomTrace : nfactors);
  Ptr<TraceSampler<Real>> traceSampler;
  if (useRandomTrace) traceSampler = makePtr<Radamacher<Real>>(theta,size);
  else                traceSampler = makePtr<TraceSampler<Real>>(theta);
  const int one(1);
  Real val = (useRandomTrace ? one/static_cast<Real>(nRandomTrace) : one);
  std::vector<Real> weight(size,val);
  bool useDeletion = list.sublist("OED").sublist("Greedy").get("Use Deletion",true);
  auto obj = makePtr<GreedyObjectiveI<Real>>(cov,factorsPV,Isampler,traceSampler,weight,useDeletion);
  return GreedyAlgorithm<Real>(obj,cvec,budget,list,os);
}

// Specialized Greedy Algorithm for I-Optimality
template<typename Real>
Ptr<Vector<Real>> GreedyAlgorithmI(const Ptr<Constraint<Real>>& model,
                                   const Ptr<Vector<Real>>& obs,
                                   const Ptr<SampleGenerator<Real>>& sampler,
                                   const Ptr<Vector<Real>>& theta,
                                   const Ptr<MomentOperator<Real>>& cov,
                                   const Ptr<Objective<Real>>& Ifun,
                                   const Ptr<SampleGenerator<Real>>& Isampler,
                                   const Ptr<Vector<Real>>& cvec,
                                         Real budget,
                                         ParameterList& list,
                                         std::ostream& os=std::cout) {
  list.sublist("OED").set("Optimality Type","I");
  cov->generateFactors(model,theta,obs,sampler);
  auto factorsPV = makePtr<Factors<Real>>(Ifun,theta,Isampler);
  bool useRandomTrace = list.sublist("OED").sublist("I-Optimality").get("Randomized Trace Estimation",false);
  int nRandomTrace    = list.sublist("OED").sublist("I-Optimality").get("Number of Samples",100);
  int nfactors        = theta->dimension();
  int size            = (useRandomTrace ? nRandomTrace : nfactors);
  Ptr<TraceSampler<Real>> traceSampler;
  if (useRandomTrace) traceSampler = makePtr<Radamacher<Real>>(theta,size);
  else                traceSampler = makePtr<TraceSampler<Real>>(theta);
  const int one(1);
  Real val = (useRandomTrace ? one/static_cast<Real>(nRandomTrace) : one);
  std::vector<Real> weight(size,val);
  bool useDeletion = list.sublist("OED").sublist("Greedy").get("Use Deletion",true);
  auto obj = makePtr<GreedyObjectiveI<Real>>(cov,factorsPV,Isampler,traceSampler,weight,useDeletion);
  return GreedyAlgorithm<Real>(obj,cvec,budget,list,os);
}

// Specialized Greedy Algorithm for I-Optimality
template<typename Real>
Ptr<Vector<Real>> GreedyAlgorithmI(const Ptr<Objective<Real>>& model,
                                   const Ptr<SampleGenerator<Real>>& sampler,
                                   const Ptr<Vector<Real>>& theta,
                                   const Ptr<MomentOperator<Real>>& cov,
                                   const Ptr<Objective<Real>>& Ifun,
                                   const Ptr<SampleGenerator<Real>>& Isampler,
                                   const Ptr<Vector<Real>>& cvec,
                                         Real budget,
                                         ParameterList& list,
                                         std::ostream& os=std::cout) {
  list.sublist("OED").set("Optimality Type","I");
  cov->generateFactors(model,theta,sampler);
  auto factorsPV = makePtr<Factors<Real>>(Ifun,theta,Isampler);
  bool useRandomTrace = list.sublist("OED").sublist("I-Optimality").get("Randomized Trace Estimation",false);
  int nRandomTrace    = list.sublist("OED").sublist("I-Optimality").get("Number of Samples",100);
  int nfactors        = theta->dimension();
  int size            = (useRandomTrace ? nRandomTrace : nfactors);
  Ptr<TraceSampler<Real>> traceSampler;
  if (useRandomTrace) traceSampler = makePtr<Radamacher<Real>>(theta,size);
  else                traceSampler = makePtr<TraceSampler<Real>>(theta);
  const int one(1);
  Real val = (useRandomTrace ? one/static_cast<Real>(nRandomTrace) : one);
  std::vector<Real> weight(size,val);
  bool useDeletion = list.sublist("OED").sublist("Greedy").get("Use Deletion",true);
  auto obj = makePtr<GreedyObjectiveI<Real>>(cov,factorsPV,Isampler,traceSampler,weight,useDeletion);
  return GreedyAlgorithm<Real>(obj,cvec,budget,list,os);
}

// Specialized Greedy Algorithm for I-Optimality
template<typename Real>
Ptr<Vector<Real>> GreedyAlgorithmI(const Ptr<Constraint<Real>>& model,
                                   const Ptr<Vector<Real>>& obs,
                                   const Ptr<SampleGenerator<Real>>& sampler,
                                   const Ptr<Vector<Real>>& theta,
                                   const Ptr<MomentOperator<Real>>& cov,
                                   const Ptr<Vector<Real>>& cvec,
                                         Real budget,
                                         ParameterList& list,
                                         std::ostream& os=std::cout) {
  return GreedyAlgorithmI<Real>(model,obs,sampler,theta,cov,model,obs,sampler,cvec,budget,list,os);
}

// Specialized Greedy Algorithm for I-Optimality
template<typename Real>
Ptr<Vector<Real>> GreedyAlgorithmI(const Ptr<Objective<Real>>& model,
                                   const Ptr<SampleGenerator<Real>>& sampler,
                                   const Ptr<Vector<Real>>& theta,
                                   const Ptr<MomentOperator<Real>>& cov,
                                   const Ptr<Vector<Real>>& cvec,
                                         Real budget,
                                         ParameterList& list,
                                         std::ostream& os=std::cout) {
  return GreedyAlgorithmI<Real>(model,sampler,theta,cov,model,sampler,cvec,budget,list,os);
}

//Ptr<ProbabilityVector<Real>> GreedyAlgorithm(const Ptr<Constraint<Real>>& model,
//                                             const Ptr<Vector<Real>>& obs,
//                                             const Ptr<SampleGenerator<Real>>& sampler,
//                                             const Ptr<Vector<Real>>& theta,
//                                             const Ptr<MomentOperator<Real>>& cov,
//                                             const Ptr<Vector<Real>>& cvec,
//                                                   Real budget,
//                                                   ParameterList& list,
//                                                   std::ostream& os=std::cout,
//                                             const Ptr<Vector<Real>>& CoptVec=nullPtr,
//                                             const Ptr<SampleGenerator<Real>>& Isampler=nullPtr,
//                                             const Ptr<Constraint<Real>>& predFun=nullPtr,
//                                             const Ptr<Vector<Real>>& predVec=nullPtr) {
//  // Set up algorithm for either adding or deleting experiments
//  bool useDeletion = list.sublist("OED").sublist("Greedy").get("Use Deletion",true);
//  Real scale = useDeletion ? static_cast<Real>(-1) : static_cast<Real>(1);
//  Real vval  = useDeletion ? static_cast<Real>(0)  : static_cast<Real>(1);
//
//  // Set up greedy objective function object
//  cov->generateFactors(model,theta,obs,sampler);
//  std::string type = list.sublist("OED").get("Optimality Type","C");
//  Ptr<GreedyObjective<Real>> obj;
//  if (type=="A") {
//    bool useRandomTrace = list.sublist("OED").sublist("A-Optimality").get("Randomized Trace Estimation",false);
//    int nRandomTrace    = list.sublist("OED").sublist("A-Optimality").get("Number of Samples",100);
//    int nfactors        = theta->dimension();
//    int size            = (useRandomTrace ? nRandomTrace : nfactors);
//    Ptr<TraceSampler<Real>> traceSampler;
//    if (useRandomTrace) traceSampler = makePtr<Radamacher<Real>>(theta,size);
//    else                traceSampler = makePtr<TraceSampler<Real>>(theta);
//    const int one(1);
//    Real val = (useRandomTrace ? one/static_cast<Real>(nRandomTrace) : one);
//    std::vector<Real> weight(size,val);
//    obj = makePtr<GreedyObjectiveA<Real>>(cov,traceSampler,weight,useDeletion);
//  }
//  else if (type=="C") {
//    Ptr<Vector<Real>> c = CoptVec;
//    if (c==nullPtr) {
//      c = theta->dual().clone();
//      c->setScalar(static_cast<Real>(1));
//    }
//    obj = makePtr<GreedyObjectiveC<Real>>(cov,c,useDeletion);
//  }
//  else if (type=="D") {
//    obj = makePtr<GreedyObjectiveD<Real>>(cov,useDeletion);
//  }
//  else if (type=="I") {
//    auto PF = predFun;
//    auto PV = predVec;
//    auto IS = Isampler;
//    if (PF == nullPtr || PV == nullPtr) {
//      PF = model;
//      PV = makePtr<SingletonVector<Real>>(1);
//    }
//    if (IS == nullPtr) IS = sampler;
//    auto factorsPV = makePtr<Factors<Real>>(PF,theta,PV,IS);
//    bool useRandomTrace = list.sublist("OED").sublist("I-Optimality").get("Randomized Trace Estimation",false);
//    int nRandomTrace    = list.sublist("OED").sublist("I-Optimality").get("Number of Samples",100);
//    int nfactors        = theta->dimension();
//    int size            = (useRandomTrace ? nRandomTrace : nfactors);
//    Ptr<TraceSampler<Real>> traceSampler;
//    if (useRandomTrace) traceSampler = makePtr<Radamacher<Real>>(theta,size);
//    else                traceSampler = makePtr<TraceSampler<Real>>(theta);
//    const int one(1);
//    Real val = (useRandomTrace ? one/static_cast<Real>(nRandomTrace) : one);
//    std::vector<Real> weight(size,val);
//    obj = makePtr<GreedyObjectiveI<Real>>(cov,factorsPV,IS,traceSampler,weight,useDeletion);
//  }
//
//  // Get parallel processing information
//  const auto& bman = sampler->getBatchManager();
//  const unsigned nbatch = sampler->numBatches();
//  const unsigned rank = sampler->batchID();
//  const unsigned mysize = sampler->numMySamples();
//  const unsigned size = sampler->numGlobalSamples();
//
//  auto pdata = makePtr<std::vector<Real>>(mysize,static_cast<Real>(1)-vval);
//  auto pvec  = makePtr<ProbabilityVector<Real>>(pdata,bman);
//  const auto& cdata = staticPtrCast<ProbabilityVector<Real>>(cvec)->getVector();
//
//  // Set up index vector
//  std::vector<unsigned> myind(size);
//  for (unsigned i=0; i<size; ++i) myind[i]=i;
//
//  // Print output to stream
//  os << std::endl;
//  if (useDeletion)
//    os << "Greedy Deletion Algorithm for Binary " << type << "-Optimality OED" << std::endl;
//  else
//    os << "Greedy Addition Algorithm for Binary " << type << "-Optimality OED" << std::endl;
//  os << "  ";
//  os << std::setw(6)  << std::left << "iter";
//  os << std::setw(15) << std::left << "value";
//  os << std::setw(15) << std::left << "cost";
//  os << std::setw(15) << std::left << "delta";
//  os << std::setw(10) << std::left << "index";
//  os << std::endl;
//
//  // Run greedy deletion algorithm
//  unsigned jind(0), i(0), mind(0);
//  Real val(0), vi(0), vij(0), cost(cvec->dot(*pvec));
//  std::vector<Real> gvi(2*nbatch), tmp(2);
//  for (i=0u; i<size; ++i) {
//    obj->update(*pvec,i);
//    val = obj->value(*pvec);
//
//    // Print output to stream
//    os << "  ";
//    os << std::setw(6)  << std::left << i;
//    os << std::setw(15) << std::left << val;
//    os << std::setw(15) << std::left << cost;
//    if (i==0u) {
//      os << std::setw(15) << std::left << "---";
//      os << std::setw(10) << std::left << "---";
//    }
//    else {
//      os << std::setw(15) << std::left << vi;
//      os << std::setw(10) << std::left << mind;
//    }
//    os << std::endl;
//
//    // Compute effectivity index
//    for (unsigned j=0u; j<myind.size(); ++j) {
//      vij = obj->computeIncrement(*pvec,myind[j]) / (*cdata)[myind[j]];
//      if (j==0u || scale*vi < scale*vij) {
//        vi = vij;
//        jind = j;
//      }
//    }
//
//    // Communicate index accross processors
//    tmp[0] = vi;
//    tmp[1] = static_cast<Real>(jind);
//    bman->gather(tmp.data(),2,gvi.data(),2*nbatch,0);
//    if (rank==0u) {
//      tmp[0] = gvi[0];
//      tmp[1] = gvi[1];
//      for (unsigned j=1u; j<nbatch; ++j) {
//        if (scale*tmp[0] < scale*gvi[2*j]) {
//          tmp[0] = gvi[2*j];
//          tmp[1] = gvi[2*j+1];
//        }
//      }
//    }
//    bman->broadcast(tmp.data(),2,0);
//    vi = tmp[0];
//    jind = static_cast<unsigned>(tmp[1]);
//    mind = myind[jind];
//
//    // Update design and cost
//    if (!useDeletion && cost+(*cdata)[mind] > budget) break;
//    (*pdata)[mind] = static_cast<Real>(vval); // Remove experiment
//    myind.erase(myind.begin()+jind);          // Remove experiment
//    cost += scale*(*cdata)[mind];             // Update cost
//    if (useDeletion && cost <= budget) break;
//  }
//
//  // Print output to stream
//  if (useDeletion) {
//    i++;
//    obj->update(*pvec,i);
//    val = obj->value(*pvec);
//    os << "  ";
//    os << std::setw(6)  << std::left << i;
//    os << std::setw(15) << std::left << val;
//    os << std::setw(15) << std::left << cost;
//    os << std::setw(15) << std::left << vi;
//    os << std::setw(10) << std::left << mind;
//    os << std::endl;
//  }
//
//  return pvec;
//}

} // End ROL::OED Namespace

#endif
