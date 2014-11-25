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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <iostream>

// Teuchos
#include <Teuchos_StandardCatchMacros.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_Factory.hpp>

#include <MueLu_RAPFactory.hpp>
#include <MueLu_SaPFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_TransPFactory.hpp>

#include <MueLu_UseDefaultTypes.hpp>

// This example demonstrates how to reuse some parts of a classical SA multigrid setup between runs.
//
// In this example, we suppose that the pattern of the matrix does not change between runs so that:
// - Aggregates can be reused
// - The tentative prolongator of Smoothed-Aggregation does not change (as it derived directly from the aggregate information).
// - The pattern of coarse grid A can be reused during its computation
//
// The resulting preconditioners are identical to multigrid preconditioners built without recycling the parts described above.
// This can be verified by using the --no-recycling option.

int main(int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  using Teuchos::TimeMonitor;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);

  bool success = false;
  bool verbose = true;
  try {
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
    Teuchos::FancyOStream& out = *fancy;

    typedef Teuchos::ScalarTraits<SC> STS;
    SC one = STS::one();

    // =========================================================================
    // Parameters initialization
    // =========================================================================
    Teuchos::CommandLineProcessor clp(false);

    GO nx = 100, ny = 100, nz = 100;
    Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D"); // manage parameters of the test case
    Xpetra::Parameters             xpetraParameters(clp);                          // manage parameters of Xpetra

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
    }

    Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();
    ParameterList galeriList = galeriParameters.GetParameterList();

    // =========================================================================
    // Problem construction
    // =========================================================================
    // For comments, see Driver.cpp
    out << "========================================================\n" << xpetraParameters << galeriParameters;
    std::string matrixType = galeriParameters.GetMatrixType();

    RCP<const Map>   map;
    RCP<MultiVector> coordinates;
    if (matrixType == "Laplace1D") {
      map = Galeri::Xpetra::CreateMap<LO, GO, Node>(lib, "Cartesian1D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("1D", map, galeriList);

    } else if (matrixType == "Laplace2D" || matrixType == "Star2D" || matrixType == "BigStar2D" || matrixType == "Elasticity2D") {
      map = Galeri::Xpetra::CreateMap<LO, GO, Node>(lib, "Cartesian2D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("2D", map, galeriList);

    } else if (matrixType == "Laplace3D" || matrixType == "Brick3D" || matrixType == "Elasticity3D") {
      map = Galeri::Xpetra::CreateMap<LO, GO, Node>(lib, "Cartesian3D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("3D", map, galeriList);
    }

    // Expand map to do multiple DOF per node for block problems
    if (matrixType == "Elasticity2D" || matrixType == "Elasticity3D")
      map = Xpetra::MapFactory<LO,GO,Node>::Build(map, (matrixType == "Elasticity2D" ? 2 : 3));

    out << "Processor subdomains in x direction: " << galeriList.get<int>("mx") << std::endl
        << "Processor subdomains in y direction: " << galeriList.get<int>("my") << std::endl
        << "Processor subdomains in z direction: " << galeriList.get<int>("mz") << std::endl
        << "========================================================" << std::endl;

    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(galeriParameters.GetMatrixType(), map, galeriList);

    RCP<Matrix> A = Pr->BuildMatrix();

    RCP<MultiVector> nullspace;
    if (matrixType == "Elasticity2D" || matrixType == "Elasticity3D") {
      nullspace = Pr->BuildNullspace();
      A->SetFixedBlockSize((matrixType == "Elasticity2D") ? 2 : 3);
    }

    // =========================================================================
    // Setups and solves
    // =========================================================================
    RCP<Vector> X = VectorFactory::Build(map);
    RCP<Vector> B = VectorFactory::Build(map);
    B->setSeed(846930886);
    B->randomize();

    const int nIts = 9;

    std::string thickSeparator = "==========================================================================================================================";
    std::string thinSeparator  = "--------------------------------------------------------------------------------------------------------------------------";

    int verbosityLevel = MueLu::Medium;

    RCP<TimeMonitor> tm;

    // =========================================================================
    // Solve #1 (no reuse)
    // =========================================================================
    out << thickSeparator << std::endl;
    {
      Hierarchy H(A);
      H.SetDefaultVerbLevel(verbosityLevel);

      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Solve #1: no reuse")));
      H.Setup();
      tm = Teuchos::null;

      X->putScalar(Teuchos::as<Scalar>(0));
      H.Iterate(*B, *X, nIts);
      out << "||Residual|| = " << Utils::ResidualNorm(*A, *X, *B)[0] << std::endl;
    }

    // =========================================================================
    // Solve #2 (reuse tentative P)
    // =========================================================================
    out << thickSeparator << std::endl;
    {
      Hierarchy H(A);
      H.SetDefaultVerbLevel(verbosityLevel);

      FactoryManager M;
      RCP<Factory> PtentFact = rcp(new TentativePFactory());
      M.SetFactory("Ptent", PtentFact);
      H.Keep      ("P",     PtentFact.get());

      // Original setup
      A->SetMaxEigenvalueEstimate(-one);
      H.Setup(M);

      H.Clear();
#ifdef HAVE_MUELU_DEBUG
      M.ResetDebugData();
#endif

      out << thinSeparator << "\nPreconditioner status:" << std::endl;
      H.print(out, MueLu::Extreme);

      // Reuse setup
      A->SetMaxEigenvalueEstimate(-one);
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Solve #2: reuse tentative P")));
      // Need same number of levels, otherwise we get in trouble when requesting the data for numLevels+1, something like this
      //  MueLu::Level(0)::GetFactory(Aggregates, 0): No FactoryManager
      //    during request for data "     Aggregates" on level 0 by factory TentativePFactory
      //    during request for data "      Nullspace" on level 1 by factory NullspaceFactory
      //    during request for data "      Nullspace" on level 1 by factory TentativePFactory
      //    during request for data "      Nullspace" on level 2 by factory NullspaceFactory
      //    during request for data "      Nullspace" on level 2 by factory TentativePFactory
      //    during request for data "              P" on level 3 by factory SaPFactory
      //    during request for data "              P" on level 3 by factory NoFactory
      H.Setup(M, 0, H.GetNumLevels());
      tm = Teuchos::null;

      X->putScalar(Teuchos::as<Scalar>(0));
      H.Iterate(*B, *X, nIts);
      out << "||Residual|| = " << Utils::ResidualNorm(*A, *X, *B)[0] << std::endl;
    }

    // =========================================================================
    // Solve #3 (reuse smoothed P and R)
    // =========================================================================
    out << thickSeparator << std::endl;
    {
      Hierarchy H(A);
      H.SetDefaultVerbLevel(verbosityLevel);

      FactoryManager M;
      RCP<Factory> PFact = rcp(new SaPFactory());
      M.SetFactory("P", PFact);
      H.Keep      ("P", PFact.get());
      RCP<Factory> RFact = rcp(new TransPFactory());
      M.SetFactory("R", RFact);
      H.Keep      ("R", RFact.get());
      RCP<Factory> RAPFact = rcp(new RAPFactory());
      if (lib == Xpetra::UseEpetra) {
        ParameterList RAPparams;
        RAPparams.set("Keep AP Pattern",  true);
        RAPparams.set("Keep RAP Pattern", true);
        RAPFact->SetParameterList(RAPparams);
      }
      M.SetFactory("A", RAPFact);

      // Original setup
      A->SetMaxEigenvalueEstimate(-one);
      H.Setup(M);

      H.Clear();
#ifdef HAVE_MUELU_DEBUG
      M.ResetDebugData();
#endif

      out << thinSeparator << "\nPreconditioner status:" << std::endl;
      H.print(out, MueLu::Extreme);

      // Reuse setup
      A->SetMaxEigenvalueEstimate(-one);
      tm = rcp(new TimeMonitor(*TimeMonitor::getNewTimer("Solve #3: reuse smoothed P and R")));
      // Need same number of levels (see comments above)
      H.Setup(M, 0, H.GetNumLevels());
      tm = Teuchos::null;

      X->putScalar(Teuchos::as<Scalar>(0));
      H.Iterate(*B, *X, nIts);
      out << "||Residual|| = " << Utils::ResidualNorm(*A, *X, *B)[0] << std::endl;
    }

    {
      const bool alwaysWriteLocal = false;
      const bool writeGlobalStats = true;
      const bool writeZeroTimers  = false;
      const bool ignoreZeroTimers = true;
      const std::string filter    = "Solve #";
      TimeMonitor::summarize(A->getRowMap()->getComm().ptr(), std::cout, alwaysWriteLocal, writeGlobalStats,
                             writeZeroTimers, Teuchos::Union, filter, ignoreZeroTimers);
    }

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
