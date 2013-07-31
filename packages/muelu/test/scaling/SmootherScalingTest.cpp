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
#include <iostream>

#include <Xpetra_MultiVectorFactory.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>
//

#include "MueLu.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"


#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_UseShortNames.hpp>


int main(int argc, char *argv[]) {
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;

  //
  // MPI initialization using Teuchos
  //

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  //
  // Parameters
  //

  Teuchos::CommandLineProcessor clp(false);

  GO nx, ny, nz;
  nx=50;
  ny=50;
  nz=50;
  Galeri::Xpetra::Parameters<GO> matrixParameters(clp, nx, ny, nz, "Laplace2D"); // manage parameters of the test case
  Xpetra::Parameters             xpetraParameters(clp);                          // manage parameters of Xpetra

  int  optSweeps  = 500;  clp.setOption("sweeps",                &optSweeps,  "smoother iterations (matvecs) to perform");
  bool optTimings = true; clp.setOption("timings", "notimings", &optTimings, "print timings to screen");

  std::string optSmoo = "chebyshev";
  clp.setOption("smoother",                &optSmoo,  "smoother type (chebyshev,sgs)");
/*
  // this doesn't work ... but it used to?!
  enum SMOOTYPES {CHEBYSHEV, SYMM_GAUSS_SEIDEL};
  const int numSmootherTypes = 2;
  const SMOOTYPES smootherValues[] = { CHEBYSHEV, SYMM_GAUSS_SEIDEL};
  const char *smootherNames[] = {"chebyshev", "sgs"};
  SMOOTYPES optSmoo = CHEBYSHEV;
  clp.setOption<SMOOTYPES>("smoother", &optSmoo,
                numSmootherTypes, smootherValues, smootherNames,
                "smoother type");
*/

  switch (clp.parse(argc, argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_ERROR:
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  if (comm->getRank() == 0) {
    std::cout << "========================================================" << std::endl
              << xpetraParameters << matrixParameters;
  }

  //
  // Construct the problem
  //

    TimeMonitor globalTimeMonitor(*TimeMonitor::getNewTimer("SmootherScalingTest: S - Global Time"));
    {

    RCP<Matrix> A;
    RCP<MultiVector> X;
    RCP<MultiVector> RHS;
    Level level;
    RCP<SmootherPrototype> smoother;
    {
      TimeMonitor tm(*TimeMonitor::getNewTimer("SmootherScalingTest: 1 - Matrix and vector creation"));

      RCP<const Map> map;

      // Retrieve matrix parameters (they may have been changed on the command line), and pass them to Galeri.
      // Galeri will attempt to create a square-as-possible distribution of subdomains di, e.g.,
      //                                 d1  d2  d3
      //                                 d4  d5  d6
      //                                 d7  d8  d9
      //                                 d10 d11 d12
      // A perfect distribution is only possible when the #processors is a perfect square.
      // This *will* result in "strip" distribution if the #processors is a prime number or if the factors are very different in
      // size. For example, np=14 will give a 7-by-2 distribution.
      // If you don't want Galeri to do this, specify mx or my on the galeriList.
      Teuchos::ParameterList pl = matrixParameters.GetParameterList();
      Teuchos::ParameterList galeriList;
      galeriList.set("nx", pl.get("nx", nx));
      galeriList.set("ny", pl.get("ny", ny));
      galeriList.set("nz", pl.get("nz", nz));

      if (matrixParameters.GetMatrixType() == "Laplace1D") {
        map = MapFactory::Build(xpetraParameters.GetLib(), matrixParameters.GetNumGlobalElements(), 0, comm);
      }
      else if (matrixParameters.GetMatrixType() == "Laplace2D" || matrixParameters.GetMatrixType() == "Star2D") {
        map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian2D", comm, galeriList);
      }
      else if (matrixParameters.GetMatrixType() == "Laplace3D") {
        map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian3D", comm, galeriList);
      }

      //FIXME
      if (comm->getRank() == 0) {
        GO mx = galeriList.get("mx", -1);
        GO my = galeriList.get("my", -1);
        std::cout << "Processor subdomains in x direction: " << mx << std::endl
                  << "Processor subdomains in y direction: " << my << std::endl
                  << "========================================================" << std::endl;
      }

      Teuchos::RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
          Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());
      A = Pr->BuildMatrix();

      // Define X, B
      X = MultiVectorFactory::Build(map, 1);
      RHS = MultiVectorFactory::Build(map, 1);

      X->setSeed(846930886);
      X->randomize();
      A->apply(*X, *RHS, Teuchos::NO_TRANS, (SC)1.0, (SC)0.0);
      Teuchos::Array<ST::magnitudeType> norms(1);
      RHS->norm2(norms);
      RHS->scale(1.0/norms[0]);

      RCP<MueLu::FactoryManagerBase> fm = rcp(new FactoryManager());
      level.SetFactoryManager(fm);
      level.SetLevelID(0);
      level.Set("A", A);
    }

    TimeMonitor tm(*TimeMonitor::getNewTimer("SmootherScalingTest: 2 - Smoother Setup"));
    {

      Teuchos::ParameterList ifpackList;
      std::string ifpackType;
      if (optSmoo == "chebyshev") {
        ifpackType = "CHEBYSHEV";
        ifpackList.set("chebyshev: degree", (LO) optSweeps);

      } else if (optSmoo == "sgs") {
        ifpackType = "RELAXATION";
        ifpackList.set("relaxation: sweeps", (LO) optSweeps);
        ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");

      } else {
        throw(MueLu::Exceptions::RuntimeError("bad smoother choice"));
      }
/*
      switch(optSmoo) {

        case SYMM_GAUSS_SEIDEL:
          ifpackType = "RELAXATION";
          ifpackList.set("relaxation: sweeps", (LO) optSweeps);
          break;

        case CHEBYSHEV:
          ifpackType = "CHEBYSHEV";
          ifpackList.set("chebyshev: degree", (LO) optSweeps);
          break;

        default:
          //should never get here
          throw(MueLu::Exceptions::RuntimeError("bad smoother choice"));
          break;
      }
*/
      smoother = rcp(new TrilinosSmoother(ifpackType, ifpackList));
      smoother->Setup(level);
    }


    RCP<Time> smootherTimer = TimeMonitor::getNewTimer("SmootherScalingTest: 3 - Smoother kernel");
    {

      TimeMonitor tm1(*smootherTimer);
    smoother->Apply(*X, *RHS);
/*
      H->IsPreconditioner(false);
      int numAMGIterations=1;
      H->Iterate(*B, numAMGIterations, *X);
*/
    }

  } //global time
  if (optTimings) {
    Teuchos::TableFormat &format = TimeMonitor::format();
    format.setPrecision(25);
    TimeMonitor::summarize();
  }

} //main
