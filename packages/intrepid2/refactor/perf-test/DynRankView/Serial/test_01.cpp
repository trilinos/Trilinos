// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
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
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

/** \file test_01.cpp
\brief  Unit test for the RealSpaceTools class.
\author Created by Kyungjoo Kim
*/

#include "Kokkos_Core.hpp"
#include <impl/Kokkos_Timer.hpp>

#include "Teuchos_CommandLineProcessor.hpp"

#include "Intrepid2_Types.hpp"
#include "test_01.hpp"

int main(int argc, char *argv[]) {

  Teuchos::CommandLineProcessor clp;
  clp.setDocString("Intrepid2::DynRankView_PerfTest01.\n");

  int nworkset = 8;
  clp.setOption("nworkset", &nworkset, "# of worksets");

  int C = 4096;
  clp.setOption("C", &C, "# of Cells in a workset");

  int P = 9;
  clp.setOption("P", &P, "# of Points in a Cell");

  int D = 3;
  clp.setOption("D", &D, "Dimension of point data (1 - 3)");

  bool verbose = true;
  clp.setOption("enable-verbose", "disable-verbose", &verbose, "Flag for verbose printing");

  clp.recogniseAllOptions(true);
  clp.throwExceptions(false);

  Teuchos::CommandLineProcessor::EParseCommandLineReturn r_parse= clp.parse( argc, argv );

  if (r_parse == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED) return 0;
  if (r_parse != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL  ) return -1;

  Kokkos::initialize();
  
  const int r_val = Intrepid2::Test::DynRankView_PerfTest01<double,Kokkos::Serial>(nworkset,
                                                                                   C,
                                                                                   P,
                                                                                   D,
                                                                                   verbose);
  
  Kokkos::finalize();
  
  return r_val;
}
