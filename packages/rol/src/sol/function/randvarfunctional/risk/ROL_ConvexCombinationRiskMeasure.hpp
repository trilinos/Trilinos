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

/** @ingroup risk_group
    \class ROL::ConvexCombinationRiskMeasure
    \brief Provides an interface for a convex combination of
           risk measures.

    This function provides the capability to produce a convex combination
    of risk measure, i.e.,
    \f[
       \mathcal{R}(X) = \sum_{k=1}^n \lambda_k \mathcal{R}_k(X)
    \f]
    where \f$\mathcal{R}_k\f$ are risk measures and \f$\lambda_k \ge 0\f$
    with \f$\lambda_1 + \ldots + \lambda_n = 1\f$.  In general,
    \f$\mathcal{R}\f$ is not law-invariant or coherent unless each
    \f$\mathcal{R}_k\f$ is.
*/

namespace ROL {

template<class Real>
class ConvexCombinationRiskMeasure : public RandVarFunctional<Real> {
private:
  typedef typename std::vector<Real>::size_type uint;

  std::vector<Real> lambda_;
  std::vector<ROL::Ptr<RandVarFunctional<Real> > > risk_;
  uint size_;
  std::vector<int> statVec_;

  Ptr<SampledScalar<Real>> values_;
  Ptr<SampledScalar<Real>> gradvecs_;
  Ptr<SampledVector<Real>> gradients_;
  Ptr<SampledVector<Real>> hessvecs_;

  using RandVarFunctional<Real>::g_;
  using RandVarFunctional<Real>::hv_;

  void initializeCCRM(void) {
    values_    = makePtr<SampledScalar<Real>>();
    gradvecs_  = makePtr<SampledScalar<Real>>();
    gradients_ = makePtr<SampledVector<Real>>();
    hessvecs_  = makePtr<SampledVector<Real>>();

    RandVarFunctional<Real>::setStorage(values_,gradients_);
    RandVarFunctional<Real>::setHessVecStorage(gradvecs_,hessvecs_);
    for (uint i = 0; i < size_; ++i) {
      risk_[i]->setStorage(values_,gradients_);
      risk_[i]->setHessVecStorage(gradvecs_,hessvecs_);
    }
  }

  void checkInputs(void) {
    uint lSize = lambda_.size(), rSize = risk_.size();
    ROL_TEST_FOR_EXCEPTION((lSize!=rSize),std::invalid_argument,
      ">>> ERROR (ROL::ConvexCombinationRiskMeasure): Convex combination parameter and risk measure arrays have different sizes!");
    Real sum(0), zero(0), one(1);
    for (uint i = 0; i < lSize; ++i) {
      ROL_TEST_FOR_EXCEPTION((lambda_[i]>one || lambda_[i]<zero), std::invalid_argument,
        ">>> ERROR (ROL::ConvexCombinationRiskMeasure): Element of convex combination parameter array out of range!");
      ROL_TEST_FOR_EXCEPTION(risk_[i] == ROL::nullPtr, std::invalid_argument,
        ">>> ERROR (ROL::ConvexCombinationRiskMeasure): Risk measure pointer is null!");
      sum += lambda_[i];
    }
    ROL_TEST_FOR_EXCEPTION((std::abs(sum-one) > std::sqrt(ROL_EPSILON<Real>())),std::invalid_argument,
      ">>> ERROR (ROL::ConvexCombinationRiskMeasure): Coefficients do not sum to one!");
    initializeCCRM();
  }

public:
  /** \brief Constructor.

      @param[in]     parlist is a parameter list specifying inputs

      parlist should contain sublists "SOL"->"Risk Measure"->"Convex Combination Risk Measure" and
      within the "Convex Combination Risk Measure" sublist should have the following parameters
      \li "Convex Combination Parameters" (greater than 0 and sum to 1)
      \li Sublists labeled 1 to n with risk measure definitions.
  */
  ConvexCombinationRiskMeasure(ROL::ParameterList &parlist)
    : RandVarFunctional<Real>(), size_(0) {
    ROL::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Convex Combination Risk Measure");
    // Get convex combination parameters
    lambda_ = ROL::getArrayFromStringParameter<Real>(list,"Convex Combination Parameters");

    size_ = lambda_.size();
    // Build risk measures
    statVec_.clear();
    risk_.clear(); risk_.resize(size_,ROL::nullPtr);
    for (uint i = 0; i < size_; ++i) {
      std::ostringstream convert;
      convert << i;
      std::string si = convert.str();
      ROL::ParameterList &ilist = list.sublist(si);
      std::string name = ilist.get<std::string>("Name");
      ROL::ParameterList riskList;
      riskList.sublist("SOL").sublist("Risk Measure").set("Name",name);
      riskList.sublist("SOL").sublist("Risk Measure").sublist(name) = ilist;
      risk_[i] = RiskMeasureFactory<Real>(riskList);
      // Get statistic information
      int nstat;
      std::vector<Real> lower, upper;
      bool isBound;
      RiskMeasureInfo(riskList,name,nstat,lower,upper,isBound);
      statVec_.push_back(nstat);
    }
    // Check inputs
    checkInputs();
  }

  void setSample(const std::vector<Real> &point, const Real weight) {
    RandVarFunctional<Real>::setSample(point,weight);
    for (uint i = 0; i < size_; ++i) {
      risk_[i]->setSample(point,weight);
    }
  }

  void resetStorage(bool flag = true) {
    RandVarFunctional<Real>::resetStorage(flag);
    for (uint i = 0; i < size_; ++i) {
      risk_[i]->resetStorage(flag);
    }
  }

  void initialize(const Vector<Real> &x) {
    RandVarFunctional<Real>::initialize(x);
    for (uint i = 0; i < size_; ++i) {
      risk_[i]->initialize(x);
    }
  }

  void updateValue(Objective<Real>         &obj,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   Real                    &tol) {
    std::vector<Real> statx;
    int offset(0);
    for (uint i = 0; i < size_; ++i) {
      statx.resize(statVec_[i]);
      for (int j = 0; j < statVec_[i]; ++j) {
        statx[j] = xstat[offset+j];
      }
      risk_[i]->updateValue(obj,x,statx,tol);
      offset += statVec_[i];
    }
  }

  Real getValue(const Vector<Real>      &x,
                const std::vector<Real> &xstat,
                SampleGenerator<Real>   &sampler) {
    Real val(0);
    std::vector<Real> statx;
    int offset(0);
    for (uint i = 0; i < size_; ++i) {
      statx.resize(statVec_[i]);
      for (int j = 0; j < statVec_[i]; ++j) {
        statx[j] = xstat[offset+j];
      }
      val += lambda_[i]*risk_[i]->getValue(x,statx,sampler);
      offset += statVec_[i];
    }
    return val;
  }

  void updateGradient(Objective<Real>         &obj,
                      const Vector<Real>      &x,
                      const std::vector<Real> &xstat,
                      Real                    &tol) {
    std::vector<Real> statx;
    int offset(0);
    for (uint i = 0; i < size_; ++i) {
      statx.resize(statVec_[i]);
      for (int j = 0; j < statVec_[i]; ++j) {
        statx[j] = xstat[offset+j];
      }
      risk_[i]->updateGradient(obj,x,statx,tol);
      offset += statVec_[i];
    }
  }

  void getGradient(Vector<Real>            &g,
                   std::vector<Real>       &gstat,
                   const Vector<Real>      &x,
                   const std::vector<Real> &xstat,
                   SampleGenerator<Real>   &sampler) {
    std::vector<Real> statg, statx;
    int offset(0);
    for (uint i = 0; i < size_; ++i) {
      statg.resize(statVec_[i]);
      statx.resize(statVec_[i]);
      for (int j = 0; j < statVec_[i]; ++j) {
        statg[j] = static_cast<Real>(0);
        statx[j] = xstat[offset+j];
      }
      g_->zero();
      risk_[i]->getGradient(*g_,statg,x,statx,sampler);
      g.axpy(lambda_[i],*g_);
      for (int j = 0; j < statVec_[i]; ++j) {
        gstat[offset+j] = lambda_[i]*statg[j];
      }
      offset += statVec_[i];
    }
  }

  void updateHessVec(Objective<Real>         &obj,
                     const Vector<Real>      &v,
                     const std::vector<Real> &vstat,
                     const Vector<Real>      &x,
                     const std::vector<Real> &xstat,
                     Real                    &tol) {
    std::vector<Real> statx, statv;
    int offset(0);
    for (uint i = 0; i < size_; ++i) {
      statx.resize(statVec_[i]);
      statv.resize(statVec_[i]);
      for (int j = 0; j < statVec_[i]; ++j) {
        statx[j] = xstat[offset+j];
        statv[j] = vstat[offset+j];
      }
      risk_[i]->updateHessVec(obj,v,statv,x,statx,tol);
      offset += statVec_[i];
    }
  }

  void getHessVec(Vector<Real>            &hv,
                  std::vector<Real>       &hvstat,
                  const Vector<Real>      &v,
                  const std::vector<Real> &vstat,
                  const Vector<Real>      &x,
                  const std::vector<Real> &xstat,
                  SampleGenerator<Real>   &sampler) {
    std::vector<Real> stath, statx, statv;
    int offset(0);
    for (uint i = 0; i < size_; ++i) {
      stath.resize(statVec_[i]);
      statx.resize(statVec_[i]);
      statv.resize(statVec_[i]);
      for (int j = 0; j < statVec_[i]; ++j) {
        stath[j] = static_cast<Real>(0);
        statx[j] = xstat[offset+j];
        statv[j] = vstat[offset+j];
      }
      hv_->zero();
      risk_[i]->getHessVec(*hv_,stath,v,statv,x,statx,sampler);
      hv.axpy(lambda_[i],*hv_);
      for (int j = 0; j < statVec_[i]; ++j) {
        hvstat[offset+j] = lambda_[i]*stath[j];
      }
      offset += statVec_[i];
    }
  }
};

}

#endif
