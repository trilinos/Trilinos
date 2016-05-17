
#ifndef BUILDSAMPLERHPP
#define BUILDSAMPLERHPP

#include "Teuchos_Comm.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// SOL Inputs
#include "ROL_MonteCarloGenerator.hpp"
#include "ROL_DistributionFactory.hpp"
#include "ROL_TpetraTeuchosBatchManager.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

template<class Real>
class BuildSampler {
private:
  Teuchos::RCP<ROL::BatchManager<Real> > bman_;
  Teuchos::RCP<ROL::SampleGenerator<Real> > sampler_;
  std::vector<Teuchos::RCP<ROL::Distribution<Real> > > distVec_;

public:
  virtual ~BuildSampler() {};
  BuildSampler(const Teuchos::RCP<const Teuchos::Comm<int> > &comm,
                     Teuchos::ParameterList &Slist,
                     Teuchos::ParameterList &Elist) {
    /*** Build stochastic functionality. ***/
    int sdim = 6;
    distVec_.clear(); distVec_.resize(sdim,Teuchos::null);
    // Build distribution for volumetric force
    bool vStoch = Slist.sublist("Problem").sublist("Volumetric Force").get<bool>("Stochastic");
    if ( vStoch ) {
      Teuchos::ParameterList &magList
        = Slist.sublist("Problem").sublist("Volumetric Force").sublist("Magnitude");
      distVec_[0] = ROL::DistributionFactory<Real>(magList);
      Elist.sublist("Elasticity").set("Bodyforce Magnitude",distVec_[0]->moment(1));
      Teuchos::ParameterList &angList
        = Slist.sublist("Problem").sublist("Volumetric Force").sublist("Angle");
      distVec_[1] = ROL::DistributionFactory<Real>(angList);
      Elist.sublist("Elasticity").set("Bodyforce Angle",distVec_[1]->moment(1));
    }
    else {
      Real magCtr = Elist.sublist("Elasticity").get<Real>("Bodyforce Magnitude");
      Teuchos::ParameterList magList;
      magList.sublist("Distribution").set("Name","Dirac");
      magList.sublist("Distribution").sublist("Dirac").set("Location",magCtr);
      distVec_[0] = ROL::DistributionFactory<Real>(magList);
      Real angCtr = Elist.sublist("Elasticity").get<Real>("Bodyforce Angle");
      Teuchos::ParameterList angList;
      angList.sublist("Distribution").set("Name","Dirac");
      angList.sublist("Distribution").sublist("Dirac").set("Location",angCtr);
      distVec_[1] = ROL::DistributionFactory<Real>(angList);
    }
    // Build distribution for traction force
    bool tStoch = Slist.sublist("Problem").sublist("Traction Force").get<bool>("Stochastic");
    if ( tStoch ) {
      Teuchos::ParameterList &magList
        = Slist.sublist("Problem").sublist("Traction Force").sublist("Magnitude");
      distVec_[2] = ROL::DistributionFactory<Real>(magList);
      Elist.sublist("Elasticity").set("Traction Magnitude",distVec_[2]->moment(1));
      Teuchos::ParameterList &angList
        = Slist.sublist("Problem").sublist("Traction Force").sublist("Angle");
      distVec_[3] = ROL::DistributionFactory<Real>(angList);
      Elist.sublist("Elasticity").set("Traction Angle",distVec_[3]->moment(1));
    }
    else {
      Real magCtr = Elist.sublist("Elasticity").get<Real>("Traction Magnitude");
      Teuchos::ParameterList magList;
      magList.sublist("Distribution").set("Name","Dirac");
      magList.sublist("Distribution").sublist("Dirac").set("Location",magCtr);
      distVec_[2] = ROL::DistributionFactory<Real>(magList);
      Real angCtr = Elist.sublist("Elasticity").get<Real>("Traction Angle");
      Teuchos::ParameterList angList;
      angList.sublist("Distribution").set("Name","Dirac");
      angList.sublist("Distribution").sublist("Dirac").set("Location",angCtr);
      distVec_[3] = ROL::DistributionFactory<Real>(angList);
    }
    // Build distribution for point force
    bool pStoch = Slist.sublist("Problem").sublist("Point Force").get<bool>("Stochastic");
    if ( pStoch ) {
      Teuchos::ParameterList &magList
        = Slist.sublist("Problem").sublist("Point Force").sublist("Magnitude");
      distVec_[4] = ROL::DistributionFactory<Real>(magList);
      Elist.sublist("Elasticity").set("Pointload Magnitude",distVec_[4]->moment(1));
      Teuchos::ParameterList &angList
        = Slist.sublist("Problem").sublist("Point Force").sublist("Angle");
      distVec_[5] = ROL::DistributionFactory<Real>(angList);
      Elist.sublist("Elasticity").set("Pointload Angle",distVec_[5]->moment(1));
    }
    else {
      Real magCtr = Elist.sublist("Elasticity").get<Real>("Pointload Magnitude");
      Teuchos::ParameterList magList;
      magList.sublist("Distribution").set("Name","Dirac");
      magList.sublist("Distribution").sublist("Dirac").set("Location",magCtr);
      distVec_[4] = ROL::DistributionFactory<Real>(magList);
      Real angCtr = Elist.sublist("Elasticity").get<Real>("Pointload Angle");
      Teuchos::ParameterList angList;
      angList.sublist("Distribution").set("Name","Dirac");
      angList.sublist("Distribution").sublist("Dirac").set("Location",angCtr);
      distVec_[5] = ROL::DistributionFactory<Real>(angList);
    }
    // Build sampler
    int nsamp = Slist.sublist("Problem").get("Number of Monte Carlo Samples",100);
    if ( !vStoch && !tStoch && !pStoch ) {
      nsamp = 1;
    }
    bman_ = Teuchos::rcp(new ROL::TpetraTeuchosBatchManager<Real>(comm));
    sampler_  = Teuchos::rcp(new ROL::MonteCarloGenerator<Real>(nsamp,distVec_,bman_,false,false,0));
  }

  Teuchos::RCP<ROL::SampleGenerator<Real> >& get(void) {
    return sampler_;
  }

  Teuchos::RCP<ROL::BatchManager<Real> >& getBatchManager(void) {
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
        file << std::setw(width) << std::left << pt[0]
             << std::setw(width) << std::left << pt[1]
             << std::setw(width) << std::left << pt[2]
             << std::setw(width) << std::left << pt[3]
             << std::setw(width) << std::left << pt[4]
             << std::setw(width) << std::left << pt[5]
             << std::setw(width) << std::left << wt
             << std::endl; 
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
