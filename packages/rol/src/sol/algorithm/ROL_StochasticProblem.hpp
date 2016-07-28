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
  Teuchos::RCP<Teuchos::ParameterList> parlist_;

  Teuchos::RCP<ParametrizedObjective<Real> > ORIGINAL_obj_;
  Teuchos::RCP<Vector<Real> >                ORIGINAL_vec_;
  Teuchos::RCP<BoundConstraint<Real> >       ORIGINAL_bnd_;

  Teuchos::RCP<Objective<Real> >       obj_;
  Teuchos::RCP<Vector<Real> >          vec_;
  Teuchos::RCP<BoundConstraint<Real> > bnd_;

  Teuchos::RCP<SampleGenerator<Real> > vsampler_;
  Teuchos::RCP<SampleGenerator<Real> > gsampler_;
  Teuchos::RCP<SampleGenerator<Real> > hsampler_;

  bool setVector_;

public:
  StochasticProblem(void)
    : OptimizationProblem<Real>(),
      parlist_(Teuchos::null),
      ORIGINAL_obj_(Teuchos::null), ORIGINAL_vec_(Teuchos::null), ORIGINAL_bnd_(Teuchos::null),
      obj_(Teuchos::null), vec_(Teuchos::null), bnd_(Teuchos::null),
      vsampler_(Teuchos::null), gsampler_(Teuchos::null), hsampler_(Teuchos::null),
      setVector_(false) {}

  StochasticProblem(Teuchos::ParameterList &parlist)
    : OptimizationProblem<Real>(),
      ORIGINAL_obj_(Teuchos::null), ORIGINAL_vec_(Teuchos::null), ORIGINAL_bnd_(Teuchos::null),
      obj_(Teuchos::null), vec_(Teuchos::null), bnd_(Teuchos::null),
      vsampler_(Teuchos::null), gsampler_(Teuchos::null), hsampler_(Teuchos::null),
      setVector_(false) {
    parlist_ = Teuchos::rcpFromRef(parlist);
  }

  StochasticProblem(Teuchos::ParameterList &parlist,
                    const Teuchos::RCP<ParametrizedObjective<Real> > &obj,
                    const Teuchos::RCP<SampleGenerator<Real> > &sampler,
                    const Teuchos::RCP<Vector<Real> > &vec)
    : OptimizationProblem<Real>(),
      ORIGINAL_obj_(obj), ORIGINAL_vec_(vec), ORIGINAL_bnd_(Teuchos::null),
      obj_(Teuchos::null), vec_(Teuchos::null), bnd_(Teuchos::null),
      vsampler_(sampler), gsampler_(sampler), hsampler_(sampler),
      setVector_(false) {
    parlist_ = Teuchos::rcpFromRef(parlist);
    setObjective(obj);
    setSolutionVector(vec);
    setBoundConstraint(Teuchos::null);
  }

  StochasticProblem(Teuchos::ParameterList &parlist,
                    const Teuchos::RCP<ParametrizedObjective<Real> > &obj,
                    const Teuchos::RCP<SampleGenerator<Real> > &vsampler,
                    const Teuchos::RCP<SampleGenerator<Real> > &gsampler,
                    const Teuchos::RCP<Vector<Real> > &vec)
    : OptimizationProblem<Real>(),
      ORIGINAL_obj_(obj), ORIGINAL_vec_(vec), ORIGINAL_bnd_(Teuchos::null),
      obj_(Teuchos::null), vec_(Teuchos::null), bnd_(Teuchos::null),
      vsampler_(vsampler), gsampler_(gsampler), hsampler_(gsampler),
      setVector_(false) {
    parlist_ = Teuchos::rcpFromRef(parlist);
    setObjective(obj);
    setSolutionVector(vec);
    setBoundConstraint(Teuchos::null);
  }

  StochasticProblem(Teuchos::ParameterList &parlist,
                    const Teuchos::RCP<ParametrizedObjective<Real> > &obj,
                    const Teuchos::RCP<SampleGenerator<Real> > &vsampler,
                    const Teuchos::RCP<SampleGenerator<Real> > &gsampler,
                    const Teuchos::RCP<SampleGenerator<Real> > &hsampler,
                    const Teuchos::RCP<Vector<Real> > &vec)
    : OptimizationProblem<Real>(),
      ORIGINAL_obj_(obj), ORIGINAL_vec_(vec), ORIGINAL_bnd_(Teuchos::null),
      obj_(Teuchos::null), vec_(Teuchos::null), bnd_(Teuchos::null),
      vsampler_(vsampler), gsampler_(gsampler), hsampler_(hsampler),
      setVector_(false) {
    parlist_ = Teuchos::rcpFromRef(parlist);
    setObjective(obj);
    setSolutionVector(vec);
    setBoundConstraint(Teuchos::null);
  }

  StochasticProblem(Teuchos::ParameterList &parlist,
                    const Teuchos::RCP<ParametrizedObjective<Real> > &obj,
                    const Teuchos::RCP<SampleGenerator<Real> > &sampler,
                    const Teuchos::RCP<Vector<Real> > &vec,
                    const Teuchos::RCP<BoundConstraint<Real> > &bnd)
    : OptimizationProblem<Real>(),
      ORIGINAL_obj_(obj), ORIGINAL_vec_(vec), ORIGINAL_bnd_(bnd),
      obj_(Teuchos::null), vec_(Teuchos::null), bnd_(Teuchos::null),
      vsampler_(sampler), gsampler_(sampler), hsampler_(sampler),
      setVector_(false) {
    parlist_ = Teuchos::rcpFromRef(parlist);
    setObjective(obj);
    setSolutionVector(vec);
    setBoundConstraint(bnd);
  }

  StochasticProblem(Teuchos::ParameterList &parlist,
                    const Teuchos::RCP<ParametrizedObjective<Real> > &obj,
                    const Teuchos::RCP<SampleGenerator<Real> > &vsampler,
                    const Teuchos::RCP<SampleGenerator<Real> > &gsampler,
                    const Teuchos::RCP<Vector<Real> > &vec,
                    const Teuchos::RCP<BoundConstraint<Real> > &bnd)
    : OptimizationProblem<Real>(),
      ORIGINAL_obj_(obj), ORIGINAL_vec_(vec), ORIGINAL_bnd_(bnd),
      obj_(Teuchos::null), vec_(Teuchos::null), bnd_(Teuchos::null),
      vsampler_(vsampler), gsampler_(gsampler), hsampler_(gsampler),
      setVector_(false) {
    parlist_ = Teuchos::rcpFromRef(parlist);
    setObjective(obj);
    setSolutionVector(vec);
    setBoundConstraint(bnd);
  }

  StochasticProblem(Teuchos::ParameterList &parlist,
                    const Teuchos::RCP<ParametrizedObjective<Real> > &obj,
                    const Teuchos::RCP<SampleGenerator<Real> > &vsampler,
                    const Teuchos::RCP<SampleGenerator<Real> > &gsampler,
                    const Teuchos::RCP<SampleGenerator<Real> > &hsampler,
                    const Teuchos::RCP<Vector<Real> > &vec,
                    const Teuchos::RCP<BoundConstraint<Real> > &bnd)
    : OptimizationProblem<Real>(),
      ORIGINAL_obj_(obj), ORIGINAL_vec_(vec), ORIGINAL_bnd_(bnd),
      obj_(Teuchos::null), vec_(Teuchos::null), bnd_(Teuchos::null),
      vsampler_(vsampler), gsampler_(gsampler), hsampler_(hsampler),
      setVector_(false) {
    parlist_ = Teuchos::rcpFromRef(parlist);
    setObjective(obj);
    setSolutionVector(vec);
    setBoundConstraint(bnd);
  }

  void setParameterList(Teuchos::ParameterList &parlist) {
    parlist_ = Teuchos::rcpFromRef(parlist);
    if (ORIGINAL_obj_ != Teuchos::null) {
      setObjective(ORIGINAL_obj_);
    }
    if (ORIGINAL_vec_ != Teuchos::null) {
      setSolutionVector(ORIGINAL_vec_);
    }
    if (ORIGINAL_bnd_ != Teuchos::null) {
      setBoundConstraint(ORIGINAL_bnd_);
    }
  }

  void setValueSampleGenerator(const Teuchos::RCP<SampleGenerator<Real> > &vsampler) {
    vsampler_ = vsampler;
    if ( gsampler_ == Teuchos::null ) {
      gsampler_ = vsampler_;
    }
    if ( hsampler_ == Teuchos::null ) {
      hsampler_ = gsampler_;
    }
  }

  void setGradientSampleGenerator(const Teuchos::RCP<SampleGenerator<Real> > &gsampler) {
    gsampler_ = gsampler;
    if ( hsampler_ == Teuchos::null ) {
      hsampler_ = gsampler_;
    }
  }

  void setHessVecSampleGenerator(const Teuchos::RCP<SampleGenerator<Real> > &hsampler) {
    hsampler_ = hsampler;
  }

  void setObjective(const Teuchos::RCP<ParametrizedObjective<Real> > &obj) {
    if ( parlist_ == Teuchos::null ) {
     TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
       ">>> ERROR (ROL::StochasticProblem): parameter list not set!");
    }
    else {
      ORIGINAL_obj_ = obj;
      if ( vsampler_ == Teuchos::null ) {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
          ">>> ERROR (ROL::StochasticProblem): value sampler not set!");
      }
      else {
        // Determine Stochastic Optimization Type
        std::string type = parlist_->sublist("SOL").get("Stochastic Optimization Type","Risk Neutral");
        if ( type == "Risk Neutral" ) {
          bool storage = parlist_->sublist("SOL").get("Store Sampled Value and Gradient",true);
          obj_ = Teuchos::rcp(new RiskNeutralObjective<Real>(obj,vsampler_,gsampler_,hsampler_,storage));
        }
        else if ( type == "Risk Averse" ) {
          obj_ = Teuchos::rcp(new RiskAverseObjective<Real>(obj,*parlist_,vsampler_,gsampler_,hsampler_));
        }
        else if ( type == "BPOE" ) {
          Real order     = parlist_->sublist("SOL").sublist("BPOE").get("Moment Order",1.);
          Real threshold = parlist_->sublist("SOL").sublist("BPOE").get("Threshold",0.);
          obj_ = Teuchos::rcp(new BPOEObjective<Real>(obj,order,threshold,vsampler_,gsampler_,hsampler_)); 
        }
        else if ( type == "Mean Value" ) {
          std::vector<Real> mean = computeSampleMean(vsampler_);
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
    }
  }

  void setSolutionVector(const Teuchos::RCP<Vector<Real> > &vec) {
    if ( parlist_ == Teuchos::null ) {
     TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
       ">>> ERROR (ROL::StochasticProblem): parameter list not set!");
    }
    else {
      ORIGINAL_vec_ = vec;
      // Determine Stochastic Optimization Type
      std::string type = parlist_->sublist("SOL").get("Stochastic Optimization Type","Risk Neutral");
      if ( type == "Risk Neutral" || type == "Mean Value" ) {
        vec_ = vec;
      }
      else if ( type == "Risk Averse" ) {
        vec_ = Teuchos::rcp(new RiskVector<Real>(*parlist_,vec));
      }
      else if ( type == "BPOE" ) {
        std::vector<Real> stat(1,1);
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
  }

  void setSolutionStatistic(const Real stat) {
    if ( parlist_ == Teuchos::null ) {
     TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
       ">>> ERROR (ROL::StochasticProblem): parameter list not set!");
    }
    else {
      if ( setVector_ ) {
        // Determine Stochastic Optimization Type
        std::string type = parlist_->sublist("SOL").get("Stochastic Optimization Type","Risk Neutral");
        if ( type == "Risk Averse" || type == "BPOE" ) {
          RiskVector<Real> &x = Teuchos::dyn_cast<RiskVector<Real> >(*vec_);
          x.setStatistic(stat);
        }
        // Set OptimizationProblem data
        OptimizationProblem<Real>::setSolutionVector(vec_);
      }
    }
  }

  void setBoundConstraint(const Teuchos::RCP<BoundConstraint<Real> > &bnd) {
    if ( parlist_ == Teuchos::null ) {
     TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
       ">>> ERROR (ROL::StochasticProblem): parameter list not set!");
    }
    else {
      ORIGINAL_bnd_ = bnd;
      // Determine Stochastic Optimization Type
      std::string type = parlist_->sublist("SOL").get("Stochastic Optimization Type","Risk Neutral");
      if ( type == "Risk Neutral" || type == "Mean Value" ) {
        bnd_ = bnd;
      }
      else if ( type == "Risk Averse" || type == "BPOE" ) {
        bnd_ = Teuchos::rcp(new RiskBoundConstraint<Real>(*parlist_,bnd));
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                                   "Invalid stochastic optimization type" << type);
      }
      // Set OptimizationProblem data
      OptimizationProblem<Real>::setBoundConstraint(bnd_);
    }
  }

//  Real getSolutionStatistic(void) {
//    try {
//      return Teuchos::dyn_cast<const RiskVector<Real> >(
//               Teuchos::dyn_cast<const Vector<Real> >(*vec_)).getStatistic();
//    }
//    catch (std::exception &e) {
//      return 0.;
//    }
//  }

  Real getSolutionStatistic(void) {
    if ( parlist_ == Teuchos::null ) {
     TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
       ">>> ERROR (ROL::StochasticProblem): parameter list not set!");
    }
    else {
      try {
        const RiskVector<Real> x = Teuchos::dyn_cast<const RiskVector<Real> >(
          Teuchos::dyn_cast<const Vector<Real> >(*vec_));
        std::string type = parlist_->sublist("SOL").get("Stochastic Optimization Type","Risk Neutral");
        Real val(0);
        if ( type == "Risk Averse" ) {
          Teuchos::ParameterList &list
            = parlist_->sublist("SOL").sublist("Risk Measure");
          std::string risk = list.get("Name","CVaR");
          if ( risk == "Mixed-Quantile Quadrangle" ) {
            Teuchos::ParameterList &MQQlist = list.sublist("Mixed-Quantile Quadrangle");
            Teuchos::Array<Real> coeff
              = Teuchos::getArrayFromStringParameter<Real>(MQQlist,"Coefficient Array");
            for (int i = 0; i < coeff.size(); i++) {
              val += coeff[i]*x.getStatistic(i);
            }
          }
          else if ( risk == "Super Quantile Quadrangle" ) {
            SuperQuantileQuadrangle<Real> sqq(*parlist_);
            val = sqq.computeStatistic(*vec_);
          }
          else if ( risk == "Chebyshev-Kusuoka" ) {
            ChebyshevKusuoka<Real> sqq(*parlist_);
            val = static_cast<SingletonKusuoka<Real> >(sqq).computeStatistic(*vec_);
          }
          else if ( risk == "Singleton Kusuoka" ) {
            SingletonKusuoka<Real> sqq(*parlist_);
            val = sqq.computeStatistic(*vec_);
          }
          else if ( risk == "Quantile-Radius Quadrangle" ) {
            Real half(0.5);
            val = half*(x.getStatistic(0) + x.getStatistic(1));
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
        return 0;
      }
    }
  }

  std::vector<std::vector<Real> > checkObjectiveGradient( const Vector<Real> &d,
                                                          const bool printToStream = true,
                                                          std::ostream & outStream = std::cout,
                                                          const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                                                          const int order = 1 ) {
    Teuchos::RCP<Vector<Real> > dp = d.clone();
    dp->set(d);
    Real stat(5.1264386);
    Teuchos::RCP<Vector<Real> > D = createVector(dp,stat);
    return OptimizationProblem<Real>::checkObjectiveGradient(*D,printToStream,outStream,numSteps,order);
  }

  std::vector<std::vector<Real> > checkObjectiveHessVec( const Vector<Real> &v,
                                                         const bool printToStream = true,
                                                         std::ostream & outStream = std::cout,
                                                         const int numSteps = ROL_NUM_CHECKDERIV_STEPS,
                                                         const int order = 1 ) {
    Teuchos::RCP<Vector<Real> > vp = v.clone();
    vp->set(v);
    Real stat(3.223468906);
    Teuchos::RCP<Vector<Real> > V = createVector(vp,stat);
    return OptimizationProblem<Real>::checkObjectiveHessVec(*V,printToStream,outStream,numSteps,order);
  } 

private:
  std::vector<Real> computeSampleMean(Teuchos::RCP<SampleGenerator<Real> > &sampler) {
    // Compute mean value of inputs and set parameter in objective
    int dim = sampler->getMyPoint(0).size(), nsamp = sampler->numMySamples();
    std::vector<Real> loc(dim), mean(dim), pt(dim);
    Real wt(0);
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

  Teuchos::RCP<Vector<Real> > createVector(Teuchos::RCP<Vector<Real> > &vec, Real stat = 1) {
    if ( parlist_ == Teuchos::null ) {
     TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
       ">>> ERROR (ROL::StochasticProblem): parameter list not set!");
    }
    else {
      // Determine Stochastic Optimization Type
      std::string type = parlist_->sublist("SOL").get("Stochastic Optimization Type","Risk Neutral");
      if ( type == "Risk Neutral" || type == "Mean Value" ) {
        return vec;
      }
      else if ( type == "Risk Averse" ) {
        return Teuchos::rcp(new RiskVector<Real>(*parlist_,vec,stat));
      }
      else if ( type == "BPOE" ) {
        std::vector<Real> statistic(1,stat);
        return Teuchos::rcp(new RiskVector<Real>(vec,statistic,true));
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                                   "Invalid stochastic optimization type" << type);
      }
    }
  }
};
}
#endif
