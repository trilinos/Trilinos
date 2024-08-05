// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <iostream>
#include "Teuchos_GlobalMPISession.hpp"

#include "ROL_Types.hpp"
#include "ROL_ParameterListConverters.hpp"

int main(int argc, char *argv[]) {

  using namespace Teuchos;

  std::string infile  = "parameters.xml";
  std::string outfile = "tiered_parameters.xml";

  auto inlist = ROL::getParametersFromXmlFile(infile);

  ROL::ParameterList outlist;

  ROL::tierParameterList(outlist,*inlist);

  ROL::writeParameterListToXmlFile(outlist,outfile);


}
