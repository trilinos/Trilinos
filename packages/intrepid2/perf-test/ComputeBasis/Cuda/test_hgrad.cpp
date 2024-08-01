// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file test_01.cpp
\brief  Unit test for the RealSpaceTools class.
\author Created by Kyungjoo Kim
*/

#include <iomanip>

#include "Kokkos_Core.hpp"
#include <Kokkos_Timer.hpp>

#include "Teuchos_CommandLineProcessor.hpp"

#include "Intrepid2_Types.hpp"
#include "test_hgrad.hpp"

int main(int argc, char *argv[]) {

  Teuchos::CommandLineProcessor clp;
  clp.setDocString("Intrepid2::DynRankView_PerfTest01.\n");

  int nworkset = 8;
  clp.setOption("nworkset", &nworkset, "# of worksets");

  int C = 4096;
  clp.setOption("C", &C, "# of Cells in a workset");

  int order = 2;
  clp.setOption("order", &order, "cubature order");

  bool verbose = true;
  clp.setOption("enable-verbose", "disable-verbose", &verbose, "Flag for verbose printing");

  clp.recogniseAllOptions(true);
  clp.throwExceptions(false);

  Teuchos::CommandLineProcessor::EParseCommandLineReturn r_parse= clp.parse( argc, argv );

  if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
  if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL  ) return -1;

  Kokkos::initialize();

  if (verbose) 
    std::cout << "Testing datatype double\n";

  const int r_val_double = Intrepid2::Test::ComputeBasis_HGRAD
    <double,Kokkos::Cuda>(nworkset,
                            C,
                            order,
                            verbose);
  return r_val_double;
}
