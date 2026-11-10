// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PYROL_STREAM
#define PYROL_STREAM

#include <iostream>
#include <fstream>
#include <ROL_Ptr.hpp>

namespace ROL {

  inline std::ofstream openOfstream(std::string filename){
    std::ofstream outdata;
    outdata.open(filename);
    return outdata;
  }

  inline void closeOfstream(std::ofstream &outdata) {
    outdata.close();
  }

  // A non-owning pointer to std::cout.
  inline ROL::Ptr<std::ostream> getCout() {
    return ROL::makePtrFromRef<std::ostream>(std::cout);
  }
}

#endif // PYROL_STREAM
