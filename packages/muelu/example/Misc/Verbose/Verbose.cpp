// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_FancyOStream.hpp>

// Kokkos
#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>

// Xpetra
#include <Xpetra_Parameters.hpp>
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsOperator.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraMatrixFactory.hpp>

// MueLu
#include "MueLu_Hierarchy.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_IfpackSmoother.hpp"
#include "MueLu_SmootherFactory.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

int main(int argc, char *argv[]) {
  using Teuchos::RCP;
  using std::cout;
  using std::endl;

  RCP<Teuchos::FancyOStream> fos = getFancyOStream(Teuchos::rcpFromRef(std::cout));


  //
  //
  //
  //  fos->setShowProcRank(true);
//   fos->setOutputToRootOnly(0);
//   Teuchos::OSTab tab(fos, 1, "TEST");
//   fos->setShowLinePrefix(true);
//   //  fos->pushLinePrefix("Lvl 1: Smoother");

//   //  fos->setTabIndentStr("my line");
//   *fos << "Hello" << std::endl << "Hello";

//   return 0;


  //
  //
  //

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  Teuchos::CommandLineProcessor clp(false);
  
  Galeri::Xpetra::Parameters<GO> matrixParameters(clp); // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);       // manage parameters of xpetra

  Teuchos::EVerbosityLevel verbLevel = Teuchos::VERB_DEFAULT;
  {
    Teuchos::EVerbosityLevel values[6] = {Teuchos::VERB_DEFAULT, Teuchos::VERB_NONE, Teuchos::VERB_LOW, Teuchos::VERB_MEDIUM, Teuchos::VERB_HIGH, Teuchos::VERB_EXTREME};
    const char* names[6] = {"DEFAULT", "NONE", "LOW", "MEDIUM", "HIGH", "EXTREME"};
    clp.setOption("verbLevel", &verbLevel, 6, values, names, "Verbose level");
    
  }
  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_ERROR:
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }
  
  matrixParameters.check();
  xpetraParameters.check();

  const RCP<const Map> map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
  RCP<Operator> A = Galeri::Xpetra::CreateCrsMatrix<SC, LO, GO, Map, CrsOperator>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());

  std::cout << std::endl << std::endl;
  //
  //
  //

  Hierarchy H;
  Level & l = *H.GetLevel();
  l.Request("A", NULL); //FIXME: remove line
  l.Set("A", A, NULL); //FIXME2: remove NULL

  //
  //
  //

  Teuchos::ParameterList smootherParamList;
  smootherParamList.set("relaxation: type", "symmetric Gauss-Seidel");
  smootherParamList.set("relaxation: sweeps", (LO) 1);
  smootherParamList.set("relaxation: damping factor", (SC) 1.0);
  
  IfpackSmoother ifpackSmoo("point relaxation stand-alone", smootherParamList);
  ifpackSmoo.setObjectLabel("My Ifpack Smoother");

  SmootherFactory smooFact(rcpFromRef(ifpackSmoo), rcpFromRef(ifpackSmoo));
  //SmootherFactory smooFact(rcpFromRef(ifpackSmoo));

  //  smooFact.describe(*fos, verbLevel);

  //ifpackSmoo.describe(*fos, verbLevel);

  std::cout << std::endl << std::endl;  
  ifpackSmoo.Setup(l);
  std::cout << std::endl << std::endl;
  
  //ifpackSmoo.describe(*fos, verbLevel);
  ifpackSmoo.print(*fos, verbLevel);

  //
  //
  //
  
  return EXIT_SUCCESS;

}

// setVerbLevel(Teuchos::VERB_HIGH);

//   std::cout << "---------- description() ---------- " << std::endl
//        << ifpackSmoo.description()
//        << std::endl;

//   std::cout << "---------- << operator ---------- " << std::endl
//        << ifpackSmoo;

//   std::cout << "---------- describe(NONE) ---------- " << std::endl;
//   ifpackSmoo.describe(*fos, Teuchos::VERB_NONE);

//   std::cout << "---------- describe(LOW) ---------- " << std::endl;
//   ifpackSmoo.describe(*fos, Teuchos::VERB_LOW);

//   std::cout << "---------- describe(MEDIUM) ---------- " << std::endl;
//   ifpackSmoo.describe(*fos, Teuchos::VERB_MEDIUM);

//   std::cout << "---------- describe(HIGH) ---------- " << std::endl;
//   ifpackSmoo.describe(*fos, Teuchos::VERB_HIGH);

//   std::cout << "---------- describe(EXTREME) ---------- " << std::endl;
//   ifpackSmoo.describe(*fos, Teuchos::VERB_EXTREME);

//   std::cout << "-----------" << std::endl;
