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

#ifndef ROL_CONVEXCOMBINATIONRISKMEASURE_HPP
#define ROL_CONVEXCOMBINATIONRISKMEASURE_HPP

#include "ROL_RiskMeasureFactory.hpp"

namespace ROL {

template<class Real>
class ConvexCombinationRiskMeasure : public RiskMeasure<Real> {
private:
  std::vector<Real> lambda_;
  std::vector<Teuchos::ParameterList> parlist_;
  std::vector<Teuchos::RCP<RiskMeasure<Real> > > risk_;
  int size_;

  Teuchos::RCP<Vector<Real> > dualVector0_;
  bool firstReset_;

  void checkInputs(void) const {
    int lSize = lambda_.size(), rSize = risk_.size();
    TEUCHOS_TEST_FOR_EXCEPTION((lSize!=rSize),std::invalid_argument,
      ">>> ERROR (ROL::ConvexCombinationRiskMeasure): Convex combination parameter and risk measure arrays have different sizes!");
    Real sum(0), zero(0), one(1);
    for (int i = 0; i < lSize; ++i) {
      TEUCHOS_TEST_FOR_EXCEPTION((lambda_[i]>one || lambda_[i]<zero), std::invalid_argument,
        ">>> ERROR (ROL::ConvexCombinationRiskMeasure): Element of convex combination parameter array out of range!");
      TEUCHOS_TEST_FOR_EXCEPTION(risk_[i] == Teuchos::null, std::invalid_argument,
        ">>> ERROR (ROL::ConvexCombinationRiskMeasure): Risk measure pointer is null!");
      sum += lambda_[i];
    }
    TEUCHOS_TEST_FOR_EXCEPTION((std::abs(sum-one) > std::sqrt(ROL_EPSILON<Real>())),std::invalid_argument,
      ">>> ERROR (ROL::ConvexCombinationRiskMeasure): Coefficients do not sum to one!");
  }

public:
  ConvexCombinationRiskMeasure(Teuchos::ParameterList &parlist)
    : RiskMeasure<Real>(), size_(0), firstReset_(true) {
    Teuchos::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Convex Combination Risk Measure");
    // Get convex combination parameters
    Teuchos::Array<Real> lambda
      = Teuchos::getArrayFromStringParameter<Real>(list,"Convex Combination Parameters");
    lambda_ = lambda.toVector();
    size_ = lambda_.size();
    // Build risk measures
    risk_.clear(); risk_.resize(size_,Teuchos::null);
    parlist_.clear(); parlist_.resize(size_);
    for (int i = 0; i < size_; ++i) {
      std::ostringstream convert;
      convert << i;
      std::string si = convert.str();
      Teuchos::ParameterList &ilist = list.sublist(si);
      std::string name = ilist.get<std::string>("Name");
      parlist_[i].sublist("SOL").sublist("Risk Measure").set("Name",name);
      parlist_[i].sublist("SOL").sublist("Risk Measure").sublist(name) = ilist;
      risk_[i] = RiskMeasureFactory<Real>(parlist_[i]);
    }
    // Check inputs
    checkInputs();
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x) {
    std::vector<Real> stat, stati;
    int N = 0, Ni = 0;
    // Must make x a risk vector with appropriate statistic
    const RiskVector<Real> &xr = Teuchos::dyn_cast<const RiskVector<Real> >(x);
    Teuchos::RCP<const Vector<Real> > xptr = xr.getVector();
    xr.getStatistic(stat);
    x0 = Teuchos::rcp_const_cast<Vector<Real> >(xptr);
    for (int i = 0; i < size_; ++i) {
      // Build temporary risk vector
      RiskVector<Real> xri(parlist_[i],x0);
      // Set statistic from original risk vector
      xri.getStatistic(stati);
      Ni = stati.size();
      for (int j = 0; j < Ni; ++j) {
        stati[j] = stat[N+j];
      }
      xri.setStatistic(stati);
      N += Ni;
      // Reset current risk measure
      risk_[i]->reset(x0,xri);
    }
    if (firstReset_) {
      dualVector0_ = x0->dual().clone();
      firstReset_ = false;
    }
    dualVector0_->zero();
  }

  void reset(Teuchos::RCP<Vector<Real> > &x0, const Vector<Real> &x,
             Teuchos::RCP<Vector<Real> > &v0, const Vector<Real> &v) {
    ConvexCombinationRiskMeasure<Real>::reset(x0,x);
    std::vector<Real> xstat, xstati, vstat, vstati;
    int N = 0, Ni = 0;
    // Must make x and v risk vectors with appropriate statistics
    const RiskVector<Real> &xr = Teuchos::dyn_cast<const RiskVector<Real> >(x);
    const RiskVector<Real> &vr = Teuchos::dyn_cast<const RiskVector<Real> >(v);
    Teuchos::RCP<const Vector<Real> > xptr = xr.getVector();
    Teuchos::RCP<const Vector<Real> > vptr = vr.getVector();
    x0 = Teuchos::rcp_const_cast<Vector<Real> >(xptr);
    v0 = Teuchos::rcp_const_cast<Vector<Real> >(vptr);
    xr.getStatistic(xstat);
    vr.getStatistic(vstat);
    for (int i = 0; i < size_; ++i) {
      // Build temporary risk vector
      RiskVector<Real> xri(parlist_[i],x0), vri(parlist_[i],v0);
      // Set statistic from original risk vector
      xri.getStatistic(xstati);
      vri.getStatistic(vstati);
      Ni = xstati.size();
      for (int j = 0; j < Ni; ++j) {
        xstati[j] = xstat[N+j];
        vstati[j] = vstat[N+j];
      }
      xri.setStatistic(xstati);
      vri.setStatistic(vstati);
      N += Ni;
      // Reset current risk measure
      risk_[i]->reset(x0,xri,v0,vri);
    }
    if (firstReset_) {
      dualVector0_ = x0->dual().clone();
      firstReset_ = false;
    }
    dualVector0_->zero();
  }

  void update(const Real val, const Real weight) {
    for (int i = 0; i < size_; ++i) {
      risk_[i]->update(val,weight);
    }
  }

  Real getValue(SampleGenerator<Real> &sampler) {
    Real val(0);
    for (int i = 0; i < size_; ++i) {
      val += lambda_[i]*risk_[i]->getValue(sampler);
    }
    return val;
  }

  void update(const Real val, const Vector<Real> &g, const Real weight) {
    for (int i = 0; i < size_; ++i) {
      risk_[i]->update(val,g,weight);
    }
  }

  void getGradient(Vector<Real> &g, SampleGenerator<Real> &sampler) {
    g.zero();
    // g does not have the correct dimension if it is a risk vector
    RiskVector<Real> &gr = Teuchos::dyn_cast<RiskVector<Real> >(g);
    std::vector<Real> stat, stati;
    for (int i = 0; i < size_; ++i) {
      RiskVector<Real> gri(parlist_[i],dualVector0_);
      risk_[i]->getGradient(gri,sampler);
      (gr.getVector())->axpy(lambda_[i],*dualVector0_);
      gri.getStatistic(stati);
      for (int j = 0; j < stati.size(); ++j) {
        stat.push_back(lambda_[i]*stati[j]);
      }
    }
    gr.setStatistic(stat);
  }

  void update(const Real val, const Vector<Real> &g, const Real gv, const Vector<Real> &hv,
              const Real weight) {
    for (int i = 0; i < size_; ++i) {
      risk_[i]->update(val,g,gv,hv,weight);
    }
  }

  void getHessVec(Vector<Real> &hv, SampleGenerator<Real> &sampler) {
    hv.zero();
    // hv does not have the correct dimension if it is a risk vector
    RiskVector<Real> &hvr = Teuchos::dyn_cast<RiskVector<Real> >(hv);
    std::vector<Real> stat, stati;
    for (int i = 0; i < size_; ++i) {
      RiskVector<Real> hvri(parlist_[i],dualVector0_);
      risk_[i]->getHessVec(hvri,sampler);
      (hvr.getVector())->axpy(lambda_[i],*dualVector0_);
      hvri.getStatistic(stati);
      for (int j = 0; j < stati.size(); ++j) {
        stat.push_back(lambda_[i]*stati[j]);
      }
    }
    hvr.setStatistic(stat);
  }
};

}

#endif
