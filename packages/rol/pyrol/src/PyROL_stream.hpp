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
#include "Teuchos_RCP.hpp"

namespace ROL {

  std::ofstream openOfstream(std::string filename){
    std::ofstream outdata;
    outdata.open(filename);
    return outdata;
  }

  void closeOfstream(std::ofstream &outdata) {
    outdata.close();
  }

  Teuchos::RCP<std::ostream> getCout() {
    Teuchos::RCP<std::ostream> out = Teuchos::rcp(&std::cout, false);
    return out;
  }
}

#endif // PYROL_STREAM
