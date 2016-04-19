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

#ifndef ROL_RISKMEASUREINFO_HPP
#define ROL_RISKMEASUREINFO_HPP

#include "Teuchos_ParameterList.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template<class Real>
inline void RiskMeasureInfo(Teuchos::ParameterList &parlist, std::string &name,
                            int &nStatistic, std::vector<Real> &lower,
                            std::vector<Real> &upper, bool isBoundActivated) {
  name = parlist.sublist("SOL").sublist("Risk Measure").get("Name","CVaR");
  Real zero(0);
  lower.clear(); upper.clear();
  nStatistic = 0; isBoundActivated = false;
  if ( name == "CVaR"                           ||
       name == "HMCR"                           ||
       name == "Moreau-Yosida CVaR"             ||
       name == "Log-Exponential Quadrangle"     ||
       name == "Log-Quantile Quadrangle"        ||
       name == "Mean-Variance Quadrangle"       ||
       name == "Quantile-Based Quadrangle"      ||
       name == "Smoothed Worst-Case Quadrangle" ||
       name == "Truncated Mean Quadrangle" ) {
    nStatistic = 1;
    lower.resize(nStatistic,ROL_NINF<Real>());
    upper.resize(nStatistic,ROL_INF<Real>());
  }
  else if ( name == "Quantile-Radius Quadrangle" ) {
    nStatistic = 2;
    lower.resize(nStatistic,ROL_NINF<Real>());
    upper.resize(nStatistic,ROL_INF<Real>());
  }
  else if ( name == "Coherent Exponential Utility" ||
            name == "KL Divergence" ) {
    nStatistic = 1;
    isBoundActivated = true;
    lower.resize(nStatistic,zero);
    upper.resize(nStatistic,ROL_INF<Real>());
  }
  else if ( name == "Chi-Squared Divergence" ) {
    nStatistic = 2;
    isBoundActivated = true;
    lower.resize(nStatistic,ROL_NINF<Real>()); lower[0] = zero;
    upper.resize(nStatistic,ROL_INF<Real>());
  }
  else if ( name == "Mixed-Quantile Quadrangle" ) {
    Teuchos::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Mixed-Quantile Quadrangle");
    Teuchos::Array<Real> prob
      = Teuchos::getArrayFromStringParameter<Real>(list,"Probability Array");
    nStatistic = prob.size();
    lower.resize(nStatistic,ROL_NINF<Real>());
    upper.resize(nStatistic,ROL_INF<Real>());
  }
  else if ( name == "Exponential Utility"             ||
            name == "Mean Plus Deviation From Target" ||
            name == "Mean Plus Deviation"             ||
            name == "Mean Plus Variance From Target"  ||
            name == "Mean Plus Variance" ) {
    nStatistic = 0;
  }
  else if ( name == "Convex Combination Risk Measure" ) {
    Teuchos::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Convex Combination Risk Measure");
    // Get convex combination parameters
    Teuchos::Array<Real> lambda
      = Teuchos::getArrayFromStringParameter<Real>(list,"Convex Combination Parameters");
    // Build risk measures
    std::vector<std::string> riskString;
    for (typename Teuchos::Array<Real>::size_type i = 0; i < lambda.size(); ++i) {
      std::ostringstream convert;
      convert << i;
      std::string si = convert.str();
      Teuchos::ParameterList &ilist = list.sublist(si);
      std::string name = ilist.get<std::string>("Name");
      riskString.push_back(name);
    }
    for (typename std::vector<Real>::size_type i = 0; i < riskString.size(); ++i) {
      if ( riskString[i] == "CVaR"                           ||
           riskString[i] == "HMCR"                           ||
           riskString[i] == "Moreau-Yosida CVaR"             ||
           riskString[i] == "Log-Exponential Quadrangle"     ||
           riskString[i] == "Log-Quantile Quadrangle"        ||
           riskString[i] == "Mean-Variance Quadrangle"       ||
           riskString[i] == "Quantile-Based Quadrangle"      ||
           riskString[i] == "Smoothed Worst-Case Quadrangle" ||
           riskString[i] == "Truncated Mean Quadrangle" ) {
        nStatistic += 1;
        lower.push_back(ROL_NINF<Real>());
        upper.push_back(ROL_INF<Real>());
      }
      else if ( riskString[i] == "Quantile-Radius Quadrangle" ) {
        nStatistic += 2;
        lower.push_back(ROL_NINF<Real>()); lower.push_back(ROL_NINF<Real>());
        upper.push_back(ROL_INF<Real>());  upper.push_back(ROL_INF<Real>());
      }
      else if ( riskString[i] == "Coherent Exponential Utility" ||
                riskString[i] == "KL Divergence" ) {
        nStatistic += 1;
        isBoundActivated = true;
        lower.push_back(zero);
        upper.push_back(ROL_INF<Real>());
      }
      else if ( riskString[i] == "Chi-Squared Divergence" ) {
        nStatistic += 2;
        isBoundActivated = true;
        lower.push_back(zero);            lower.push_back(ROL_NINF<Real>());
        upper.push_back(ROL_INF<Real>()); upper.push_back(ROL_INF<Real>());
      }
      else if ( riskString[i] == "Mixed-Quantile Quadrangle" ) {
        Teuchos::ParameterList &MQlist = list.sublist("Mixed-Quantile Quadrangle");
        Teuchos::Array<Real> prob
          = Teuchos::getArrayFromStringParameter<Real>(MQlist,"Probability Array");
        nStatistic += prob.size();
        for (typename Teuchos::Array<Real>::size_type j = 0; j < prob.size(); ++j) {
          lower.push_back(ROL_NINF<Real>());
          upper.push_back(ROL_INF<Real>());
        }
      }
      else if ( riskString[i] == "Exponential Utility"             ||
                riskString[i] == "Mean Plus Deviation From Target" ||
                riskString[i] == "Mean Plus Deviation"             ||
                riskString[i] == "Mean Plus Variance From Target"  ||
                riskString[i] == "Mean Plus Variance" ) {
        nStatistic += 0;
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
          ">>> (ROL::RiskMeasureInfo): Invalid risk measure " << riskString[i] << "!");
      }
    }
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
      ">>> (ROL::RiskMeasureInfo): Invalid risk measure " << name << "!");
  }
}

}
#endif
