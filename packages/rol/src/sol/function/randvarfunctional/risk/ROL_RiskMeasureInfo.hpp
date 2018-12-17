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

#include "ROL_ParameterList.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template<class Real>
inline void RiskMeasureInfo(ROL::ParameterList &parlist, std::string &name,
                            int &nStatistic, std::vector<Real> &lower,
                            std::vector<Real> &upper, bool &isBoundActivated,
                            const bool printToStream = false,
                            std::ostream &outStream = std::cout) {
  name = parlist.sublist("SOL").sublist("Risk Measure").get<std::string>("Name");
  Real zero(0);
  lower.clear(); upper.clear();
  nStatistic = 0; isBoundActivated = false;
  if ( name == "CVaR"                           ||
       name == "HMCR"                           ||
       name == "Moreau-Yosida CVaR"             ||
       name == "Generalized Moreau-Yosida CVaR" ||
       name == "Log Quantile"                   ||
       name == "Smoothed Worst Case"            ||
       name == "Safety Margin"                  ||
       name == "Log Exponential"                ||
       name == "Truncated Mean" ) {
    nStatistic = 1;
    lower.resize(nStatistic,ROL_NINF<Real>());
    upper.resize(nStatistic,ROL_INF<Real>());
  }
  else if ( name == "Quantile Radius" ) {
    nStatistic = 2;
    lower.resize(nStatistic,ROL_NINF<Real>());
    upper.resize(nStatistic,ROL_INF<Real>());
  }
  else if ( name == "Coherent Entropic Risk" ||
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
  else if ( name == "Mixed CVaR" ) {
    ROL::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Mixed CVaR");
    std::vector<Real> prob
      = ROL::getArrayFromStringParameter<Real>(list,"Probability Array");
    nStatistic = prob.size();
    lower.resize(nStatistic,ROL_NINF<Real>());
    upper.resize(nStatistic,ROL_INF<Real>());
  }
  else if ( name == "Second Order CVaR"       ||
            name == "Chebyshev Spectral Risk" ||
            name == "Spectral Risk" ) {
    ROL::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist(name);
    nStatistic = list.get("Number of Quadrature Points",5);
    lower.resize(nStatistic,ROL_NINF<Real>());
    upper.resize(nStatistic,ROL_INF<Real>());
  }
  else if ( name == "Entropic Risk"                        ||
            name == "Mean Plus Semi-Deviation From Target" ||
            name == "Mean Plus Semi-Deviation"             ||
            name == "Mean Plus Deviation From Target"      ||
            name == "Mean Plus Deviation"                  ||
            name == "Mean Plus Variance From Target"       ||
            name == "Mean Plus Variance" ) {
    nStatistic = 0;
  }
  else if ( name == "Convex Combination Risk Measure" ) {
    ROL::ParameterList &list
      = parlist.sublist("SOL").sublist("Risk Measure").sublist("Convex Combination Risk Measure");
    // Get convex combination parameters
    std::vector<Real> lambda
      = ROL::getArrayFromStringParameter<Real>(list,"Convex Combination Parameters");
    // Build risk measures
    std::vector<std::string> riskString;
    for (typename std::vector<Real>::size_type i = 0; i < lambda.size(); ++i) {
      std::ostringstream convert;
      convert << i;
      std::string si = convert.str();
      ROL::ParameterList &ilist = list.sublist(si);
      std::string namei = ilist.get<std::string>("Name");
      riskString.push_back(namei);
    }
    for (typename std::vector<Real>::size_type i = 0; i < riskString.size(); ++i) {
      if ( riskString[i] == "CVaR"                           ||
           riskString[i] == "HMCR"                           ||
           riskString[i] == "Moreau-Yosida CVaR"             ||
           riskString[i] == "Generalized Moreau-Yosida CVaR" ||
           riskString[i] == "Log Quantile"                   ||
           riskString[i] == "Smoothed Worst Case"            ||
           riskString[i] == "Safety Margin"                  ||
           riskString[i] == "Log Exponential"                ||
           riskString[i] == "Truncated Mean" ) {
        nStatistic += 1;
        lower.push_back(ROL_NINF<Real>());
        upper.push_back(ROL_INF<Real>());
      }
      else if ( riskString[i] == "Quantile Radius" ) {
        nStatistic += 2;
        lower.push_back(ROL_NINF<Real>()); lower.push_back(ROL_NINF<Real>());
        upper.push_back(ROL_INF<Real>());  upper.push_back(ROL_INF<Real>());
      }
      else if ( riskString[i] == "Coherent Entropic Risk" ||
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
      else if ( riskString[i] == "Mixed CVaR" ) {
        ROL::ParameterList &MQlist = list.sublist("Mixed CVaR");
        std::vector<Real> prob
          = ROL::getArrayFromStringParameter<Real>(MQlist,"Probability Array");
        nStatistic += prob.size();
        for (typename std::vector<Real>::size_type j = 0; j < prob.size(); ++j) {
          lower.push_back(ROL_NINF<Real>());
          upper.push_back(ROL_INF<Real>());
        }
      }
      else if ( riskString[i] == "Second Order CVaR"       ||
                riskString[i] == "Chebyshev Spectral Risk" ||
                riskString[i] == "Spectral Risk" ) {
        ROL::ParameterList &SQlist = list.sublist(riskString[i]);
        int nSQQstat = SQlist.get("Number of Quadrature Points",5);
        nStatistic += nSQQstat;
        for (int j = 0; j < nSQQstat; ++j) {
          lower.push_back(ROL_NINF<Real>());
          upper.push_back(ROL_INF<Real>());
        }
      }
      else if ( riskString[i] == "Entropic Risk"                        ||
                riskString[i] == "Mean Plus Semi-Deviation From Target" ||
                riskString[i] == "Mean Plus Semi-Deviation"             ||
                riskString[i] == "Mean Plus Deviation From Target"      ||
                riskString[i] == "Mean Plus Deviation"                  ||
                riskString[i] == "Mean Plus Variance From Target"       ||
                riskString[i] == "Mean Plus Variance" ) {
        nStatistic += 0;
      }
      else {
        ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
          ">>> (ROL::RiskMeasureInfo): Invalid risk measure " << riskString[i] << "!");
      }
    }
  }
  else {
    ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
      ">>> (ROL::RiskMeasureInfo): Invalid risk measure " << name << "!");
  }

  // Print Information
  if ( printToStream ) {
    ROL::nullstream oldFormatState;
    oldFormatState.copyfmt(outStream);

    outStream << std::endl;
    outStream << std::scientific << std::setprecision(6);
    outStream << std::setfill('-') << std::setw(80) << "-" << std::endl;
    outStream << "  RISK MEASURE INFORMATION" << std::endl;
    outStream << std::setfill('-') << std::setw(80) << "-" << std::endl;
    outStream << "  NAME" << std::endl;
    outStream << "    " << name << std::endl;
    outStream << "  NUMBER OF STATISTICS" << std::endl;
    outStream << "    " << nStatistic << std::endl;
    outStream << "  ARE BOUNDS ACTIVATED" << std::endl;
    outStream << "    " << (isBoundActivated ? "TRUE" : "FALSE") << std::endl;
    if ( isBoundActivated ) {
      outStream << "  STATISTIC LOWER BOUNDS" << std::endl;
      for (int i = 0; i < nStatistic-1; ++i) {
        outStream << "    " << lower[i] << std::endl;
      }
      outStream << "    " << lower[nStatistic-1] << std::endl;
      outStream << "  STATISTIC UPPER BOUNDS" << std::endl;
      for (int i = 0; i < nStatistic-1; ++i) {
        outStream << "    " << upper[i] << std::endl;
      }
      outStream << "    " << upper[nStatistic-1] << std::endl;
    }
    outStream << std::setfill('-') << std::setw(80) << "-" << std::endl;
    outStream << std::endl;

    outStream.copyfmt(oldFormatState);
  }
}

}
#endif
