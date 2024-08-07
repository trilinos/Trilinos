// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BUILDSAMPLERHPP
#define BUILDSAMPLERHPP

#include "Teuchos_Comm.hpp"
#include "ROL_Stream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// SOL Inputs
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_SROMGenerator.hpp"
#include "ROL_DistributionFactory.hpp"
#include "ROL_TpetraTeuchosBatchManager.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

template<class Real>
class BuildSampler {
private:
  ROL::Ptr<ROL::BatchManager<Real> > bman_;
  ROL::Ptr<ROL::SampleGenerator<Real> > sampler_;
  std::vector<ROL::Ptr<ROL::Distribution<Real> > > distVec_;

public:
  virtual ~BuildSampler() {};
  BuildSampler(const ROL::Ptr<const Teuchos::Comm<int> > &comm,
                     Teuchos::ParameterList &Slist,
                     Teuchos::ParameterList &Elist) {
    /*** Build stochastic functionality. ***/
    Real scale = static_cast<Real>(M_PI)/static_cast<Real>(180);
    int sdim = 0;
    distVec_.clear();
    // Build distribution for volumetric force
    bool vStochMag = Elist.sublist("Elasticity").get<bool>("Stochastic Bodyforce Magnitude");
    if ( vStochMag ) {
      Teuchos::ParameterList &magList
        = Slist.sublist("Problem").sublist("Volumetric Force").sublist("Magnitude");
      distVec_.push_back(ROL::DistributionFactory<Real>(magList));
      Elist.sublist("Elasticity").set("Bodyforce Magnitude",distVec_[sdim]->moment(1));
      sdim++;
    }
    bool vStochAng = Elist.sublist("Elasticity").get<bool>("Stochastic Bodyforce Angle");
    if ( vStochAng ) {
      Teuchos::ParameterList &angList
        = Slist.sublist("Problem").sublist("Volumetric Force").sublist("Angle");
      distVec_.push_back(ROL::DistributionFactory<Real>(angList));
      Elist.sublist("Elasticity").set("Bodyforce Angle",scale * distVec_[sdim]->moment(1));
      sdim++;
    }
    // Build distribution for traction force
    bool tStochMag = Elist.sublist("Elasticity").get<bool>("Stochastic Traction Magnitude");
    if ( tStochMag ) {
      Teuchos::ParameterList &magList
        = Slist.sublist("Problem").sublist("Traction Force").sublist("Magnitude");
      distVec_.push_back(ROL::DistributionFactory<Real>(magList));
      Elist.sublist("Elasticity").set("Traction Magnitude",distVec_[sdim]->moment(1));
      sdim++;
    }
    bool tStochAng = Elist.sublist("Elasticity").get<bool>("Stochastic Traction Angle");
    if ( tStochAng ) {
      Teuchos::ParameterList &angList
        = Slist.sublist("Problem").sublist("Traction Force").sublist("Angle");
      distVec_.push_back(ROL::DistributionFactory<Real>(angList));
      Elist.sublist("Elasticity").set("Traction Angle",scale * distVec_[sdim]->moment(1));
      sdim++;
    }
    // Build distribution for point force
    bool pStochMag = Elist.sublist("Elasticity").get<bool>("Stochastic Pointload Magnitude");
    if ( pStochMag ) {
      Teuchos::ParameterList &magList
        = Slist.sublist("Problem").sublist("Point Force").sublist("Magnitude");
      distVec_.push_back(ROL::DistributionFactory<Real>(magList));
      Elist.sublist("Elasticity").set("Pointload Magnitude",distVec_[sdim]->moment(1));
      sdim++;
    }
    bool pStochAng = Elist.sublist("Elasticity").get<bool>("Stochastic Pointload Angle");
    if ( pStochAng ) {
      Teuchos::ParameterList &angList
        = Slist.sublist("Problem").sublist("Point Force").sublist("Angle");
      distVec_.push_back(ROL::DistributionFactory<Real>(angList));
      Elist.sublist("Elasticity").set("Pointload Angle",scale * distVec_[sdim]->moment(1));
      sdim++;
    }
    // Build sampler
    int nsamp = Slist.sublist("Problem").get("Number of Samples",100);
    if ( !vStochMag && !vStochAng &&
         !tStochMag && !tStochAng && 
         !pStochMag && !pStochAng ) {
      nsamp = 1;
    }
    bman_ = ROL::makePtr<ROL::TpetraTeuchosBatchManager<Real>>(comm);
    bool useOBS = Slist.sublist("Problem").get("Use Optimization-Based Sampling",false);
    if ( useOBS ) {
      Slist.sublist("SOL").sublist("Sample Generator").sublist("SROM").set("Number of Samples",nsamp);
      sampler_ = ROL::makePtr<ROL::SROMGenerator<Real>>(Slist,bman_,distVec_);
    }
    else {
      sampler_ = ROL::makePtr<ROL::MonteCarloGenerator<Real>>(nsamp,distVec_,bman_,false,false,0);
    }
  }

  ROL::Ptr<ROL::SampleGenerator<Real> >& get(void) {
    return sampler_;
  }

  ROL::Ptr<ROL::BatchManager<Real> >& getBatchManager(void) {
    return bman_;
  }

  void print(const std::string &filename, const int prec = 12) const {
    int width = prec + 5 + 4;
    std::stringstream name;
    name << filename << "_" << sampler_->batchID() << ".txt";
    std::ofstream file(name.str().c_str());
    if (file.is_open()) {
      file << std::scientific << std::setprecision(prec);
      for (int i = 0; i < sampler_->numMySamples(); ++i) {
        std::vector<Real> pt = sampler_->getMyPoint(i);
        Real wt = sampler_->getMyWeight(i);
        for (int j = 0; j < static_cast<int>(pt.size()); ++j) {
          file << std::setw(width) << std::left << pt[j];
        }
        file << std::setw(width) << std::left << wt << std::endl; 
      }
      file.close();
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        ">>> (BuildSampler::print): Unable to open file!");
    }
  }

};

#endif
