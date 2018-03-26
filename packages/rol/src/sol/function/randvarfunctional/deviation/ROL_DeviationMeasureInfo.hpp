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

#ifndef ROL_DEVIATIONMEASUREINFO_HPP
#define ROL_DEVIATIONMEASUREINFO_HPP

#include "Teuchos_ParameterList.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template<class Real>
inline void DeviationMeasureInfo(Teuchos::ParameterList &parlist, std::string &name,
                            int &nStatistic, std::vector<Real> &lower,
                            std::vector<Real> &upper, bool &isBoundActivated,
                            const bool printToStream = false,
                            std::ostream &outStream = std::cout) {
  name = parlist.sublist("SOL").sublist("Deviation Measure").get<std::string>("Name");
  lower.clear(); upper.clear();
  nStatistic = 0; isBoundActivated = false;
  if ( name == "CVaR"                           ||
       name == "Moreau-Yosida CVaR"             ||
       name == "Generalized Moreau-Yosida CVaR" ||
       name == "Variance"                       ||
       name == "Entropic"                       ||
       name == "Log Quantile"                   ||
       name == "Smoothed Upper Range"           ||
       name == "Truncated Mean" ) {
    nStatistic = 1;
    lower.resize(nStatistic,ROL_NINF<Real>());
    upper.resize(nStatistic,ROL_INF<Real>());
  }
  else {
    TEUCHOS_TEST_FOR_EXCEPTION(true,std::invalid_argument,
      ">>> (ROL::DeviationMeasureInfo): Invalid deviation measure " << name << "!");
  }

  // Print Information
  if ( printToStream ) {
    Teuchos::oblackholestream oldFormatState;
    oldFormatState.copyfmt(outStream);

    outStream << std::endl;
    outStream << std::scientific << std::setprecision(6);
    outStream << std::setfill('-') << std::setw(80) << "-" << std::endl;
    outStream << "  DEVIATION MEASURE INFORMATION" << std::endl;
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
