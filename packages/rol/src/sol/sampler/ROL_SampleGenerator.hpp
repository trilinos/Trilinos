// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_SAMPLEGENERATOR_HPP
#define ROL_SAMPLEGENERATOR_HPP

#include "ROL_BatchManager.hpp"
#include "ROL_Vector.hpp"
#include <fstream>

namespace ROL {

template<class Real>
class SampleGenerator {
private:
  int begin_;
  ROL::Ptr<BatchManager<Real> > bman_;
  std::vector<std::vector<Real> > points_;
  std::vector<Real> weights_;

protected:
  void setPoints(std::vector<std::vector<Real> > &p) {
    points_.clear();
    points_.assign(p.begin(),p.end());
  }
  void setWeights(std::vector<Real> &w) {
    weights_.clear();
    weights_.assign(w.begin(),w.end());
  }

public:
  virtual ~SampleGenerator() {}
  SampleGenerator(const ROL::Ptr<BatchManager<Real> > &bman)
    : begin_(0), bman_(bman) {}
  SampleGenerator(const SampleGenerator<Real> &sampler)
    : begin_(sampler.begin_), bman_(sampler.bman_),
      points_(sampler.points_), weights_(sampler.weights_) {}

  virtual void update(const Vector<Real> &x) {
    begin_ = 0;
  }

  virtual int start(void) {
    return begin_;
  }

  virtual Real computeError(std::vector<Real> &vals) {
    return 0.0;
  }

  virtual Real computeError(std::vector<ROL::Ptr<Vector<Real> > > &vals, const Vector<Real> &x) {
    return 0.0;
  }

  virtual void refine(void) {
    begin_ = numMySamples();
  }

  virtual void setSamples(bool inConstructor = false) {}

  virtual int numGlobalSamples(void) const {
    Real nms  = static_cast<Real>(numMySamples());
    Real gnms = static_cast<Real>(0);
    sumAll(&nms,&gnms,1);
    return static_cast<int>(gnms);
  }

  virtual int numMySamples(void) const {
    return weights_.size();
  }

  virtual std::vector<Real> getMyPoint(const int i) const {
    return points_[i];
  }

  virtual Real getMyWeight(const int i) const {
    return weights_[i];
  }

  int batchID(void) const {
    return bman_->batchID();
  }

  int numBatches(void) const {
    return bman_->numBatches();
  }

  void sumAll(Real *input, Real *output, int dim) const {
    bman_->sumAll(input, output, dim);
  }

  void sumAll(Vector<Real> &input, Vector<Real> &output) const {
    bman_->sumAll(input,output);
  }

  void broadcast(Real *input, int cnt, int root) const {
    bman_->broadcast(input,cnt,root);
  }

  void barrier(void) const {
    bman_->barrier();
  }

  const ROL::Ptr<BatchManager<Real>> getBatchManager(void) const {
    return bman_;
  }

  void print(const std::string &filename = "samples",
             const int prec = 12) const {
    int width = prec + 5 + 4;
    std::stringstream name;
    name << filename << "_" << batchID() << ".txt";
    std::ofstream file(name.str().c_str());
    if (file.is_open()) {
      file << std::scientific << std::setprecision(prec);
      for (int i = 0; i < numMySamples(); ++i) {
        std::vector<Real> pt = getMyPoint(i);
        Real wt = getMyWeight(i);
        for (int j = 0; j < static_cast<int>(pt.size()); ++j) {
          file << std::setw(width) << std::left << pt[j];
        }
        file << std::setw(width) << std::left << wt << std::endl;
      }
      file.close();
    }
    else {
      ROL_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> (ROL::SampleGenerator::print): Unable to open file!");
    }
  }

};

}

#endif
