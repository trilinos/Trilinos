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

#ifndef ROL_STOCHASTICPROBLEM_HPP
#define ROL_STOCHASTICPROBLEM_HPP

#include "ROL_OptimizationProblem.hpp"
#include "ROL_SampleGenerator.hpp"

// Risk-Neutral Includes
#include "ROL_RiskNeutralObjective.hpp"

// Risk-Averse Includes
#include "ROL_RiskAverseObjective.hpp"
#include "ROL_RiskVector.hpp"
#include "ROL_RiskBoundConstraint.hpp"

// BPOE Includes
#include "ROL_BPOEObjective.hpp"

#include "Teuchos_ParameterList.hpp"

namespace ROL {

template<class Real>
class StochasticProblem : public OptimizationProblem<Real> {
private:
  Teuchos::RCP<Objective<Real> >       obj_;
  Teuchos::RCP<Vector<Real> >          vec_;
  Teuchos::RCP<BoundConstraint<Real> > con_;

  bool setVector_;

public:
  StochasticProblem(void)
    : obj_(Teuchos::null), vec_(Teuchos::null), con_(Teuchos::null),
      setVector_(false) {}

  StochasticProblem(Teuchos::ParameterList &parlist,
                    Teuchos::RCP<ParametrizedObjective<Real> > &obj,
                    Teuchos::RCP<SampleGenerator<Real> > &sampler,
                    Teuchos::RCP<Vector<Real> > &vec)
    : OptimizationProblem<Real>(),
      obj_(Teuchos::null), vec_(Teuchos::null), con_(Teuchos::null),
      setVector_(false) {
    setObjective(parlist,obj,sampler);
    setSolutionVector(parlist,vec);
  }

  StochasticProblem(Teuchos::ParameterList &parlist,
                    Teuchos::RCP<ParametrizedObjective<Real> > &obj,
                    Teuchos::RCP<SampleGenerator<Real> > &vsampler,
                    Teuchos::RCP<SampleGenerator<Real> > &gsampler,
                    Teuchos::RCP<Vector<Real> > &vec)
    : OptimizationProblem<Real>(),
      obj_(Teuchos::null), vec_(Teuchos::null), con_(Teuchos::null),
      setVector_(false) {
    setObjective(parlist,obj,vsampler,gsampler);
    setSolutionVector(parlist,vec);
  }

  StochasticProblem(Teuchos::ParameterList &parlist,
                    Teuchos::RCP<ParametrizedObjective<Real> > &obj,
                    Teuchos::RCP<SampleGenerator<Real> > &vsampler,
                    Teuchos::RCP<SampleGenerator<Real> > &gsampler,
                    Teuchos::RCP<SampleGenerator<Real> > &hsampler,
                    Teuchos::RCP<Vector<Real> > &vec)
    : OptimizationProblem<Real>(),
      obj_(Teuchos::null), vec_(Teuchos::null), con_(Teuchos::null),
      setVector_(false) {
    setObjective(parlist,obj,vsampler,gsampler,hsampler);
    setSolutionVector(parlist,vec);
  }

  StochasticProblem(Teuchos::ParameterList &parlist,
                    Teuchos::RCP<ParametrizedObjective<Real> > &obj,
                    Teuchos::RCP<SampleGenerator<Real> > &sampler,
                    Teuchos::RCP<Vector<Real> > &vec,
                    Teuchos::RCP<BoundConstraint<Real> > &con)
    : OptimizationProblem<Real>(),
      obj_(Teuchos::null), vec_(Teuchos::null), con_(Teuchos::null),
      setVector_(false) {
    setObjective(parlist,obj,sampler);
    setSolutionVector(parlist,vec);
    setBoundConstraint(parlist,con);
  }

  StochasticProblem(Teuchos::ParameterList &parlist,
                    Teuchos::RCP<ParametrizedObjective<Real> > &obj,
                    Teuchos::RCP<SampleGenerator<Real> > &vsampler,
                    Teuchos::RCP<SampleGenerator<Real> > &gsampler,
                    Teuchos::RCP<Vector<Real> > &vec,
                    Teuchos::RCP<BoundConstraint<Real> > &con)
    : OptimizationProblem<Real>(),
      obj_(Teuchos::null), vec_(Teuchos::null), con_(Teuchos::null),
      setVector_(false) {
    setObjective(parlist,obj,vsampler,gsampler);
    setSolutionVector(parlist,vec);
    setBoundConstraint(parlist,con);
  }

  StochasticProblem(Teuchos::ParameterList &parlist,
                    Teuchos::RCP<ParametrizedObjective<Real> > &obj,
                    Teuchos::RCP<SampleGenerator<Real> > &vsampler,
                    Teuchos::RCP<SampleGenerator<Real> > &gsampler,
                    Teuchos::RCP<SampleGenerator<Real> > &hsampler,
                    Teuchos::RCP<Vector<Real> > &vec,
                    Teuchos::RCP<BoundConstraint<Real> > &con)
    : OptimizationProblem<Real>(),
      obj_(Teuchos::null), vec_(Teuchos::null), con_(Teuchos::null),
      setVector_(false) {
    setObjective(parlist,obj,vsampler,gsampler,hsampler);
    setSolutionVector(parlist,vec);
    setBoundConstraint(parlist,con);
  }

  void setObjective(Teuchos::ParameterList &parlist,
                    Teuchos::RCP<ParametrizedObjective<Real> > &obj,
                    Teuchos::RCP<SampleGenerator<Real> > &vsampler) {
    // Determine Stochastic Optimization Type
    std::string type = parlist.sublist("SOL").get("Stochastic Optimization Type","Risk Neutral");
    if ( type == "Risk Neutral" ) {
      bool storage = parlist.sublist("SOL").get("Store Sampled Value and Gradient",true);
      obj_ = Teuchos::rcp(new RiskNeutralObjective<Real>(obj,vsampler,storage));
    }
    else if ( type == "Risk Averse" ) {
      obj_ = Teuchos::rcp(new RiskAverseObjective<Real>(obj,parlist,vsampler));
    }
    else if ( type == "BPOE" ) {
      Real order     = parlist.sublist("SOL").sublist("BPOE").get("Moment Order",1.);
      Real threshold = parlist.sublist("SOL").sublist("BPOE").get("Threshold",0.);
      obj_ = Teuchos::rcp(new BPOEObjective<Real>(obj,order,threshold,vsampler)); 
    }
    else if ( type == "Mean Value" ) {
      std::vector<Real> mean = computeSampleMean(vsampler);
      obj->setParameter(mean);
      obj_ = obj;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                                 "Invalid stochastic optimization type" << type);
    }
    // Set OptimizationProblem data
    OptimizationProblem<Real>::setObjective(obj_);
  }

  void setObjective(Teuchos::ParameterList &parlist,
                    Teuchos::RCP<ParametrizedObjective<Real> > &obj,
                    Teuchos::RCP<SampleGenerator<Real> > &vsampler,
                    Teuchos::RCP<SampleGenerator<Real> > &gsampler) {
    // Determine Stochastic Optimization Type
    std::string type = parlist.sublist("SOL").get("Stochastic Optimization Type","Risk Neutral");
    if ( type == "Risk Neutral" ) {
      bool storage = parlist.sublist("SOL").get("Store Sampled Value and Gradient",true);
      obj_ = Teuchos::rcp(new RiskNeutralObjective<Real>(obj,vsampler,gsampler,storage));
    }
    else if ( type == "Risk Averse" ) {
      obj_ = Teuchos::rcp(new RiskAverseObjective<Real>(obj,parlist,vsampler,gsampler));
    }
    else if ( type == "BPOE" ) {
      Real order     = parlist.sublist("SOL").sublist("BPOE").get("Moment Order",1.);
      Real threshold = parlist.sublist("SOL").sublist("BPOE").get("Threshold",0.);
      obj_ = Teuchos::rcp(new BPOEObjective<Real>(obj,order,threshold,vsampler,gsampler)); 
    }
    else if ( type == "Mean Value" ) {
      std::vector<Real> mean = computeSampleMean(vsampler);
      obj->setParameter(mean);
      obj_ = obj;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                                 "Invalid stochastic optimization type" << type);
    }
    // Set OptimizationProblem data
    OptimizationProblem<Real>::setObjective(obj_);
  }

  void setObjective(Teuchos::ParameterList &parlist,
                    Teuchos::RCP<ParametrizedObjective<Real> > &obj,
                    Teuchos::RCP<SampleGenerator<Real> > &vsampler,
                    Teuchos::RCP<SampleGenerator<Real> > &gsampler,
                    Teuchos::RCP<SampleGenerator<Real> > &hsampler) {
    // Determine Stochastic Optimization Type
    std::string type = parlist.sublist("SOL").get("Stochastic Optimization Type","Risk Neutral");
    if ( type == "Risk Neutral" ) {
      bool storage = parlist.sublist("SOL").get("Store Sampled Value and Gradient",true);
      obj_ = Teuchos::rcp(new RiskNeutralObjective<Real>(obj,vsampler,gsampler,hsampler,storage));
    }
    else if ( type == "Risk Averse" ) {
      obj_ = Teuchos::rcp(new RiskAverseObjective<Real>(obj,parlist,vsampler,gsampler,hsampler));
    }
    else if ( type == "BPOE" ) {
      Real order     = parlist.sublist("SOL").sublist("BPOE").get("Moment Order",1.);
      Real threshold = parlist.sublist("SOL").sublist("BPOE").get("Threshold",0.);
      obj_ = Teuchos::rcp(new BPOEObjective<Real>(obj,order,threshold,vsampler,gsampler,hsampler)); 
    }
    else if ( type == "Mean Value" ) {
      std::vector<Real> mean = computeSampleMean(vsampler);
      obj->setParameter(mean);
      obj_ = obj;
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                                 "Invalid stochastic optimization type" << type);
    }
    // Set OptimizationProblem data
    OptimizationProblem<Real>::setObjective(obj_);
  }

  void setSolutionVector(Teuchos::ParameterList &parlist,
                         Teuchos::RCP<Vector<Real> > &vec) {
    // Determine Stochastic Optimization Type
    std::string type = parlist.sublist("SOL").get("Stochastic Optimization Type","Risk Neutral");
    if ( type == "Risk Neutral" || type == "Mean Value" ) {
      vec_ = vec;
    }
    else if ( type == "Risk Averse" ) {
      vec_ = Teuchos::rcp(new RiskVector<Real>(parlist,vec));
    }
    else if ( type == "BPOE" ) {
      std::vector<Real> stat(1,1.);
      vec_ = Teuchos::rcp(new RiskVector<Real>(vec,stat,true));
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                                 "Invalid stochastic optimization type" << type);
    }
    // Set OptimizationProblem data
    OptimizationProblem<Real>::setSolutionVector(vec_);
    setVector_ = true;
  }

  void setSolutionStatistic(Teuchos::ParameterList &parlist,
                            Real stat) {
    if ( setVector_ ) {
      // Determine Stochastic Optimization Type
      std::string type = parlist.sublist("SOL").get("Stochastic Optimization Type","Risk Neutral");
      if ( type == "Risk Averse" || type == "BPOE" ) {
        RiskVector<Real> &x = Teuchos::dyn_cast<RiskVector<Real> >(*vec_);
        x.setStatistic(stat);
      }
      // Set OptimizationProblem data
      OptimizationProblem<Real>::setSolutionVector(vec_);
    }
  }

  void setBoundConstraint(Teuchos::ParameterList &parlist,
                          Teuchos::RCP<BoundConstraint<Real> > &con) {
    // Determine Stochastic Optimization Type
    std::string type = parlist.sublist("SOL").get("Stochastic Optimization Type","Risk Neutral");
    if ( type == "Risk Neutral" || type == "Mean Value" ) {
      con_ = con;
    }
    else if ( type == "Risk Averse" ) {
      std::string name = parlist.sublist("SOL").sublist("Risk Measure").get("Name","CVaR");
      if ( name == "KL Divergence" ) {
        con_ = Teuchos::rcp(new RiskBoundConstraint<Real>("BPOE",con));
      }
      else {
        con_ = Teuchos::rcp(new RiskBoundConstraint<Real>(parlist,con));
      }
    }
    else if ( type == "BPOE" ) {
      con_ = Teuchos::rcp(new RiskBoundConstraint<Real>("BPOE",con));
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                                 "Invalid stochastic optimization type" << type);
    }
    // Set OptimizationProblem data
    OptimizationProblem<Real>::setBoundConstraint(con_);
  }

  Teuchos::RCP<Vector<Real> > createVector(Teuchos::ParameterList &parlist,
                                           Teuchos::RCP<Vector<Real> > &vec) {
    // Determine Stochastic Optimization Type
    std::string type = parlist.sublist("SOL").get("Stochastic Optimization Type","Risk Neutral");
    if ( type == "Risk Neutral" || type == "Mean Value" ) {
      return vec;
    }
    else if ( type == "Risk Averse" ) {
      return Teuchos::rcp(new RiskVector<Real>(parlist,vec));
    }
    else if ( type == "BPOE" ) {
      std::vector<Real> stat(1,1.);
      return Teuchos::rcp(new RiskVector<Real>(vec,stat,true));
    }
    else {
      TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                                 "Invalid stochastic optimization type" << type);
    }
  }

  Real getSolutionStatistic(void) {
    try {
      return Teuchos::dyn_cast<const RiskVector<Real> >(
               Teuchos::dyn_cast<const Vector<Real> >(*vec_)).getStatistic();
    }
    catch (std::exception &e) {
      return 0.;
    }
  }

  Real getSolutionStatistic(Teuchos::ParameterList &parlist) {
    try {
      const RiskVector<Real> x = Teuchos::dyn_cast<const RiskVector<Real> >(
        Teuchos::dyn_cast<const Vector<Real> >(*vec_));
      std::string type = parlist.sublist("SOL").get("Stochastic Optimization Type","Risk Neutral");
      Real val = 0.0;
      if ( type == "Risk Averse" ) {
        Teuchos::ParameterList &list
          = parlist.sublist("SOL").sublist("Risk Measure");
        std::string risk = list.get("Name","CVaR");
        if ( risk == "Mixed-Quantile Quadrangle" ) {
          Teuchos::ParameterList &MQQlist = list.sublist("Mixed-Quantile Quadrangle");
          Teuchos::Array<Real> coeff
            = Teuchos::getArrayFromStringParameter<Real>(MQQlist,"Coefficient Array");
          for (int i = 0; i < coeff.size(); i++) {
            val += coeff[i]*x.getStatistic(i);
          }
        }
        else if ( risk == "Quantile-Radius Quadrangle" ) {
          val = 0.5*(x.getStatistic(0) + x.getStatistic(1));
        }
        else {
          val = x.getStatistic();
        }
      }
      else {
        val = x.getStatistic();
      }
      return val;
    }
    catch (std::exception &e) {
      return 0.;
    }
  }

private:
  std::vector<Real> computeSampleMean(Teuchos::RCP<SampleGenerator<Real> > &sampler) {
    // Compute mean value of inputs and set parameter in objective
    int dim = sampler->getMyPoint(0).size(), nsamp = sampler->numMySamples();
    std::vector<Real> loc(dim,0.0), mean(dim,0.0), pt(dim,0.0);
    Real wt = 0.0;
    for (int i = 0; i < nsamp; i++) {
      pt = sampler->getMyPoint(i);
      wt = sampler->getMyWeight(i);
      for (int j = 0; j < dim; j++) {
        loc[j] += wt*pt[j];
      }
    }
    sampler->sumAll(&loc[0],&mean[0],dim);
    return mean;
  }
};
}
#endif
