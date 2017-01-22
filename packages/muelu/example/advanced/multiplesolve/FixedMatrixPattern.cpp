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

#include <Teuchos_StandardCatchMacros.hpp>

#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_Factory.hpp>

#include <MueLu_CoupledAggregationFactory.hpp>
#include <MueLu_TentativePFactory.hpp>
#include <MueLu_RAPFactory.hpp>
#include <MueLu_NoFactory.hpp>

// This example demonstrates how to reuse some parts of a multigrid setup between runs.
//
// In this example, we suppose that the pattern of the matrix does not change between runs so that:
// - Aggregates can be reused
// - The tentative prolongator of Smoothed-Aggregation does not change (as it derived directly from the aggregate information).
// - The pattern of coarse grid A can be reused during its computation
//
// The resulting preconditioners are identical to multigrid preconditioners built without recycling the parts described above.
// This can be verified by using the --no-recycling option.
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP;

  typedef Teuchos::ScalarTraits<SC> ST;

  bool success = false;
  bool verbose = true;
  try {
    RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    //
    // Parameters
    //

    Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 8748);
    Xpetra::Parameters xpetraParameters(clp);

    bool optRecycling           = true;  clp.setOption("recycling",             "no-recycling",             &optRecycling,           "Enable recycling of the multigrid preconditioner");

    /* DO NOT WORK YET
       bool optRecyclingRAPpattern = true;  clp.setOption("recycling-rap-pattern", "no-recycling-rap-pattern", &optRecyclingRAPpattern, "Enable recycling of Ac=RAP pattern");
       bool optRecyclingAPpattern  = false; clp.setOption("recycling-ap-pattern",  "no-recycling-ap-pattern",  &optRecyclingAPpattern,  "Enable recycling of AP pattern");
       */
    bool optRecyclingRAPpattern = false;
    bool optRecyclingAPpattern  = false;

    switch (clp.parse(argc, argv)) {
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
      case Teuchos::CommandLineProcessor::PARSE_ERROR:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
    }

    // option dependencies
    if (optRecycling == false) {
      optRecyclingRAPpattern = false;
      optRecyclingAPpattern  = false;
    }

    //
    // Construct the problems
    //

    RCP<const Map> map = MapFactory::Build(lib, matrixParameters.GetNumGlobalElements(), 0, comm);
    Teuchos::RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(matrixParameters.GetMatrixType(), map, matrixParameters.GetParameterList());
    RCP<Matrix> A1 = Pr->BuildMatrix();

    RCP<Matrix> A2 = Pr->BuildMatrix(); // TODO: generate another problem would be more meaningful (ex: scale A1)

    //
    // First solve
    //

    FactoryManager M;

    Hierarchy H(A1);

    if (optRecycling) {
      // Configuring "Keep" options before running the first setup.

      // Note: "Keep" flags should not be set on the default factories provided by the FactoryManager because in the current
      // implementation of FactoryManager, those factories are freed and reallocated between Hierarchy::Setup() calls (see FactoryManager::Clean() and Hierarchy::Setup()).
      // So we define our own factories here.

      // AGGREGATES:
      // Note: aggregates are only used to build Ptent, so it is not useful to keep them. Keeping Ptent is enough.
      // RCP<Factory> AggFact   = rcp(new CoupledAggregationFactory());
      // M.SetFactory("Aggregates", AggFact);
      // H.Keep("Aggregates", AggFact.get());

      // PTENT:
      RCP<Factory> PtentFact = rcp(new TentativePFactory());
      M.SetFactory("Ptent", PtentFact);
      H.Keep("P",           PtentFact.get());
    }

    RCP<Factory> AcFact = rcp(new RAPFactory());
    M.SetFactory("A", AcFact);

    if (optRecyclingRAPpattern) {
      H.Keep("RAP graph", AcFact.get());
    }
    if (optRecyclingAPpattern) {
      H.Keep("AP graph", AcFact.get());
    }
    //

    H.Setup(M);

    {
      RCP<Vector> X = VectorFactory::Build(map);
      RCP<Vector> B = VectorFactory::Build(map);

      X->putScalar((Scalar) 0.0);
      B->setSeed(846930886); B->randomize();

      int nIts = 9;
      H.Iterate(*B, *X, nIts);

      typename ST::magnitudeType residualNorms = Utilities::ResidualNorm(*A1, *X, *B)[0];
      if (comm->getRank() == 0)
        std::cout << "||Residual|| = " << residualNorms << std::endl;
    }

    //
    // Second solve
    //

    std::cout << "Status of the preconditioner between runs:" << std::endl;
    H.print(*getFancyOStream(Teuchos::rcpFromRef(std::cout)), MueLu::High);

    // Change the problem
    RCP<Level> finestLevel = H.GetLevel(0);
    finestLevel->Set("A", A2);

    if (optRecycling) {
      // Optional: this makes sure that the aggregates are never requested (and built) during the second run.
      // If someone request the aggregates, an exception will be thrown.
      M.SetFactory("Aggregates", MueLu::NoFactory::getRCP());
    }

    // Redo the setup
    H.Setup(M);

    {
      RCP<Vector> X = VectorFactory::Build(map);
      RCP<Vector> B = VectorFactory::Build(map);

      X->putScalar((Scalar) 0.0);
      B->setSeed(846930886); B->randomize();

      int nIts = 9;
      H.Iterate(*B, *X, nIts);

      typename ST::magnitudeType residualNorms = Utilities::ResidualNorm(*A2, *X, *B)[0];
      if (comm->getRank() == 0)
        std::cout << "||Residual|| = " << residualNorms << std::endl;
    }

    //
    // Clean-up
    //

    // Remove kept data from the preconditioner. This will force recomputation on future runs. "Keep" flags are also removed.

    if (optRecycling) {
      //if aggregates explicitly kept: H.Delete("Aggregates", M.GetFactory("Aggregates").get());
      H.Delete("P",           M.GetFactory("Ptent").get());
    }
    if (optRecyclingRAPpattern) {
      H.Delete("RAP graph", M.GetFactory("A").get());
    }
    if (optRecyclingAPpattern) {
      H.Delete("AP graph", M.GetFactory("A").get());
    }

    std::cout << "Status of the preconditioner at the end:" << std::endl;
    H.print(*getFancyOStream(Teuchos::rcpFromRef(std::cout)), MueLu::High);

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}

int main(int argc, char* argv[]) {
  bool success = false;
  bool verbose = true;

  Teuchos::GlobalMPISession mpiSession(&argc,&argv);

  try {
    const bool throwExceptions     = false;
    const bool recogniseAllOptions = false;

    Teuchos::CommandLineProcessor clp(throwExceptions, recogniseAllOptions);
    Xpetra::Parameters xpetraParameters(clp);

    std::string node = "";  clp.setOption("node", &node, "node type (serial | openmp | cuda)");

    switch (clp.parse(argc, argv, NULL)) {
      case Teuchos::CommandLineProcessor::PARSE_ERROR:               return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
    }

    Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

    if (lib == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA
      return main_<double,int,int,Xpetra::EpetraNode>(clp, lib, argc, argv);
#else
      throw MueLu::Exceptions::RuntimeError("Epetra is not available");
#endif
    }

    if (lib == Xpetra::UseTpetra) {
#ifdef HAVE_MUELU_TPETRA
      if (node == "") {
        typedef KokkosClassic::DefaultNode::DefaultNodeType Node;

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double,int,long,Node>(clp, lib, argc, argv);
#else
#    if   defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_INT)
        return main_<double,int,int,Node> (clp, lib, argc, argv);
#  elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_LONG)
        return main_<double,int,long,Node>(clp, lib, argc, argv);
#  elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        return main_<double,int,long long,Node>(clp, lib, argc, argv);
#  else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#  endif
#endif
      } else if (node == "serial") {
#ifdef KOKKOS_HAVE_SERIAL
        typedef Kokkos::Compat::KokkosSerialWrapperNode Node;

#  ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double,int,long,Node>(clp, lib, argc, argv);
#  else
#    if   defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_INT)
        return main_<double,int,int,Node> (clp, lib, argc, argv);
#    elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG)
        return main_<double,int,long,Node>(clp, lib, argc, argv);
#    elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_SERIAL) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        return main_<double,int,long long,Node>(clp, lib, argc, argv);
#    else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#    endif
#  endif
#else
        throw MueLu::Exceptions::RuntimeError("Serial node type is disabled");
#endif
      } else if (node == "openmp") {
#ifdef KOKKOS_HAVE_OPENMP
        typedef Kokkos::Compat::KokkosOpenMPWrapperNode Node;

#  ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double,int,long,Node>(clp, lib, argc, argv);
#  else
#    if   defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_INT)
        return main_<double,int,int,Node> (clp, lib, argc, argv);
#    elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG)
        return main_<double,int,long,Node>(clp, lib, argc, argv);
#    elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_OPENMP) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        return main_<double,int,long long,Node>(clp, lib, argc, argv);
#    else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#    endif
#  endif
#else
        throw MueLu::Exceptions::RuntimeError("OpenMP node type is disabled");
#endif
      } else if (node == "cuda") {
#ifdef KOKKOS_HAVE_CUDA
        typedef Kokkos::Compat::KokkosCudaWrapperNode Node;

#  ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
        return main_<double,int,long,Node>(clp, lib, argc, argv);
#  else
#    if defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_TPETRA_INST_INT_INT)
        return main_<double,int,int,Node> (clp, lib, argc, argv);
#    elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_TPETRA_INST_INT_LONG)
        return main_<double,int,long,Node>(clp, lib, argc, argv);
#    elif defined(HAVE_TPETRA_INST_DOUBLE) && defined(HAVE_TPETRA_INST_CUDA) && defined(HAVE_TPETRA_INST_INT_LONG_LONG)
        return main_<double,int,long long,Node>(clp, lib, argc, argv);
#    else
        throw MueLu::Exceptions::RuntimeError("Found no suitable instantiation");
#    endif
#  endif
#else
        throw MueLu::Exceptions::RuntimeError("CUDA node type is disabled");
#endif
      } else {
        throw MueLu::Exceptions::RuntimeError("Unrecognized node type");
      }
#else
      throw MueLu::Exceptions::RuntimeError("Tpetra is not available");
#endif
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
