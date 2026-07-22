// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_REGRETINFO_HPP
#define ROL_REGRETINFO_HPP

#include "ROL_ParameterList.hpp"
#include "ROL_Types.hpp"

namespace ROL {

template<class Real>
inline void RegretMeasureInfo(ROL::ParameterList &parlist, std::string &name,
                            int &nStatistic, std::vector<Real> &lower,
                            std::vector<Real> &upper, bool &isBoundActivated,
                            const bool printToStream = false,
                            std::ostream &outStream = std::cout) {
  name = parlist.sublist("SOL").sublist("Regret Measure").get<std::string>("Name");
  lower.clear(); upper.clear();
  nStatistic = 0; isBoundActivated = false;
  if ( name == "Mean Absolute Loss"                           ||
       name == "Moreau-Yosida Mean Absolute Loss"             ||
       name == "Generalized Moreau-Yosida Mean Absolute Loss" ||
       name == "Exponential"                                  ||
       name == "Mean L2"                                      ||
       name == "Truncated Mean"                               ||
       name == "Log Quantile"                                 ||
       name == "Smoothed Worst Case" ) {
    nStatistic = 0;
  }
  else {
    ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
      ">>> (ROL::RegretMeasureInfo): Invalid probability " << name << "!");
  }

  // Print Information
  if ( printToStream ) {
    ROL::nullstream oldFormatState;
    oldFormatState.copyfmt(outStream);

    outStream << std::endl;
    outStream << std::scientific << std::setprecision(6);
    outStream << std::setfill('-') << std::setw(80) << "-" << std::endl;
    outStream << "  REGRET INFORMATION" << std::endl;
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
