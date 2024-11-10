// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef ROL_RANDVARFUNCTIONALINFO_HPP
#define ROL_RANDVARFUNCTIONALINFO_HPP

#include "ROL_RiskMeasureInfo.hpp"
#include "ROL_DeviationMeasureInfo.hpp"
#include "ROL_ErrorMeasureInfo.hpp"
#include "ROL_RegretMeasureInfo.hpp"
#include "ROL_ProbabilityInfo.hpp"

namespace ROL {

template<class Real>
inline void RandVarFunctionalInfo(ROL::ParameterList &parlist, std::string &name,
                                  int &nStatistic, std::vector<Real> &lower,
                                  std::vector<Real> &upper, bool &isBoundActivated,
                                  const bool printToStream = false,
                                  std::ostream &outStream = std::cout) {
  std::string type = parlist.sublist("SOL").get("Type","Risk Averse");
  if (type == "Risk Averse") {
    RiskMeasureInfo<Real>(parlist,name,nStatistic,lower,upper,isBoundActivated,printToStream,outStream);
  }
  else if (type == "Deviation") {
    DeviationMeasureInfo<Real>(parlist,name,nStatistic,lower,upper,isBoundActivated,printToStream,outStream);
  }
  else if (type == "Error") {
    ErrorMeasureInfo<Real>(parlist,name,nStatistic,lower,upper,isBoundActivated,printToStream,outStream);
  }
  else if (type == "Regret") {
    RegretMeasureInfo<Real>(parlist,name,nStatistic,lower,upper,isBoundActivated,printToStream,outStream);
  }
  else if (type == "Probability") {
    ProbabilityInfo<Real>(parlist,name,nStatistic,lower,upper,isBoundActivated,printToStream,outStream);
  }
  else {
    ROL_TEST_FOR_EXCEPTION(true,std::invalid_argument,
      ">>> (ROL::RandVarFunctionalInfo): Invalid random variable functional type!");
  }

  // Print Information
  if ( printToStream ) {
    ROL::nullstream oldFormatState;
    oldFormatState.copyfmt(outStream);

    outStream << std::endl;
    outStream << std::scientific << std::setprecision(6);
    outStream << std::setfill('-') << std::setw(80) << "-" << std::endl;
    outStream << "  RANDVARFUNCTIONAL INFORMATION" << std::endl;
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
