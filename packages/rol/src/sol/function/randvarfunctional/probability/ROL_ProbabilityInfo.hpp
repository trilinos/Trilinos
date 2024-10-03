// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_PROBABILITYINFO_HPP
#define ROL_PROBABILITYINFO_HPP

#include "ROL_ParameterList.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template<class Real>
inline void ProbabilityInfo(ROL::ParameterList &parlist, std::string &name,
                            int &nStatistic, std::vector<Real> &lower,
                            std::vector<Real> &upper, bool &isBoundActivated,
                            const bool printToStream = false,
                            std::ostream &outStream = std::cout) {
  name = parlist.sublist("SOL").sublist("Probability").get<std::string>("Name");
  Real zero(0);
  lower.clear(); upper.clear();
  nStatistic = 0; isBoundActivated = false;
  if ( name == "bPOE" ) {
    nStatistic = 1;
    isBoundActivated = true;
    lower.resize(nStatistic,zero);
    upper.resize(nStatistic,ROL_INF<Real>());
  }
  else if ( name == "Smoothed POE" ) {
    nStatistic = 0;
  }
  else {
    ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
      ">>> (ROL::ProbabilityInfo): Invalid probability " << name << "!");
  }

  // Print Information
  if ( printToStream ) {
    ROL::nullstream oldFormatState;
    oldFormatState.copyfmt(outStream);

    outStream << std::endl;
    outStream << std::scientific << std::setprecision(6);
    outStream << std::setfill('-') << std::setw(80) << "-" << std::endl;
    outStream << "  PROBABILITY INFORMATION" << std::endl;
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
