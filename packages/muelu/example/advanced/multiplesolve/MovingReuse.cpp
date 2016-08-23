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
#include <chrono>

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

#ifdef HAVE_MUELU_EXPLICIT_INSTANTIATION
#include <MueLu_ExplicitInstantiation.hpp>
#endif

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>      // => This header defines Belos::MueLuOp
#endif

// This example demonstrates how to reuse some parts of a classical SA multigrid setup between runs.
//
// In this example, we suppose that the pattern of the matrix does not change between runs so that:
// - Aggregates can be reused
// - The tentative prolongator of Smoothed-Aggregation does not change (as it derived directly from the aggregate information).
// - The pattern of coarse grid A can be reused during its computation
//
// The resulting preconditioners are identical to multigrid preconditioners built without recycling the parts described above.
// This can be verified by using the --no-recycling option.

#include <MueLu_CreateXpetraPreconditioner.hpp>


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Map, class Matrix>
Teuchos::RCP<Matrix> BuildMatrix(Teuchos::ParameterList& list, const Teuchos::RCP<const Map>& map) {
  typedef GlobalOrdinal GO;
  typedef LocalOrdinal  LO;
  typedef Scalar        SC;
  using  Teuchos::ArrayView;

  GO nx = list.get("nx", (GO) -1);
  GO ny = list.get("ny", (GO) -1);
  GO nz = list.get("nz", (GO) -1);
  if (nx == -1 || ny == -1 || nz == -1) {
    GO n = map->getGlobalNumElements();
    nx = (GO) Teuchos::ScalarTraits<double>::pow(n, 0.33334);
    ny = nx; nz = nx;
    TEUCHOS_TEST_FOR_EXCEPTION(nx * ny * nz != n, std::logic_error, "You need to specify nx, ny, and nz");
  }

  double  one = 1.0;
  SC  stretchx = list.get("stretchx", one);
  SC  stretchy = list.get("stretchy", one);
  SC  stretchz = list.get("stretchz", one);

  // bool keepBCs = list.get("keepBCs", false);

  SC b = -one / (stretchx*stretchx);
  SC c = -one / (stretchx*stretchx);
  SC d = -one / (stretchy*stretchy);
  SC e = -one / (stretchy*stretchy);
  SC f = -one / (stretchz*stretchz);
  SC g = -one / (stretchz*stretchz);
  SC a = -(b + c + d + e + f + g);

  LO nnz = 7;

  Teuchos::RCP<Matrix> A = Galeri::Xpetra::MatrixTraits<Map,Matrix>::Build(map, nnz);

  LO numMyElements = map->getNodeNumElements();
  GO indexBase     = map->getIndexBase();

  ArrayView<const GO> myGlobalElements = map->getNodeElementList();

  std::vector<GO> inds(nnz);
  std::vector<SC> vals(nnz);

  //    e
  //  b a c
  //    d
  // + f bottom and g top
  GO center, left, right, bottom, top, front, back;
  for (LO i = 0; i < numMyElements; ++i) {
    size_t n = 0;

    center = myGlobalElements[i] - indexBase;
    Galeri::Xpetra::GetNeighboursCartesian3d(center, nx, ny, nz, left, right, front, back, bottom, top);

    if (left   != -1) { inds[n] = left;   vals[n++] = b; }
    if (right  != -1) { inds[n] = right;  vals[n++] = c; }
    if (front  != -1) { inds[n] = front;  vals[n++] = d; }
    if (back   != -1) { inds[n] = back;   vals[n++] = e; }
    if (bottom != -1) { inds[n] = bottom; vals[n++] = f; }
    if (top    != -1) { inds[n] = top;    vals[n++] = g; }

    // diagonal (ignore Neumann for now => no update
    inds[n]   = center;
    vals[n++] = a;

    for (size_t j = 0; j < n; j++)
      inds[j] += indexBase;

    ArrayView<GO> iv(&inds[0], n);
    ArrayView<SC> av(&vals[0], n);
    A->insertGlobalValues(myGlobalElements[i], iv, av);
  }

  A->fillComplete();

  return A;
}


template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void ConstructData(const std::string& matrixType, Teuchos::ParameterList& galeriList,
                   Xpetra::UnderlyingLib lib, Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                   Teuchos::RCP<Xpetra::Matrix      <Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A,
                   Teuchos::RCP<const Xpetra::Map   <LocalOrdinal,GlobalOrdinal, Node> >&       map,
                   Teuchos::RCP<Xpetra::MultiVector <Scalar,LocalOrdinal,GlobalOrdinal,Node> >& coordinates,
                   Teuchos::RCP<Xpetra::MultiVector <Scalar,LocalOrdinal,GlobalOrdinal,Node> >& nullspace) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::TimeMonitor;

  map         = Galeri::Xpetra::CreateMap<LO, GO, Node>(lib, "Cartesian3D", comm, galeriList);
  coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("3D", map, galeriList);

  A = BuildMatrix<SC,LO,GO,Map,CrsMatrixWrap>(galeriList, map);
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, int argc, char *argv[]) {
#include <MueLu_UseShortNames.hpp>
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::TimeMonitor;
  using namespace std::chrono;

  // =========================================================================
  // MPI initialization using Teuchos
  // =========================================================================
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // =========================================================================
  // Convenient definitions
  // =========================================================================
  typedef Teuchos::ScalarTraits<SC> STS;
  SC one = STS::one(), zero = STS::zero();

  RCP<Teuchos::FancyOStream> fancy = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  Teuchos::FancyOStream& out = *fancy;
  out.setOutputToRootOnly(0);

  // =========================================================================
  // Parameters initialization
  // =========================================================================
  GO nx = 100, ny = 100, nz = 100;
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace2D"); // manage parameters of the test case
  Xpetra::Parameters             xpetraParameters(clp);                          // manage parameters of Xpetra

  std::string xmlFileName = "";     clp.setOption("xml",                &xmlFileName, "read parameters from a file");
  int         numRebuilds = 0;      clp.setOption("rebuild",            &numRebuilds, "#times to rebuild hierarchy");
  bool        useFilter   = true;   clp.setOption("filter", "nofilter", &useFilter,   "Print out only Setup times");
  bool        modify      = true;   clp.setOption("modify", "nomodify", &modify,      "Change values of the matrix used for reuse");
  std::string solveType   = "cg";   clp.setOption("solver",             &solveType,   "solve type: (none | cg | standalone)");
  double      tol         = 1e-6;   clp.setOption("tol",                &tol,         "solver convergence tolerance");

  clp.recogniseAllOptions(true);
  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
  }
  Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

  // Retrieve matrix parameters (they may have been changed on the command line)
  // [for instance, if we changed matrix type from 2D to 3D we need to update nz]
  ParameterList galeriList = galeriParameters.GetParameterList();

  // =========================================================================
  // Problem construction
  // =========================================================================
  // For comments, see Driver.cpp
  out << "========================================================\n" << xpetraParameters << galeriParameters;
  std::string matrixType = galeriParameters.GetMatrixType();

  out << "Processor subdomains in x direction: " << galeriList.get<GO>("mx") << std::endl
      << "Processor subdomains in y direction: " << galeriList.get<GO>("my") << std::endl
      << "Processor subdomains in z direction: " << galeriList.get<GO>("mz") << std::endl
      << "========================================================" << std::endl;

  // =========================================================================
  // Setups and solves
  // =========================================================================
  std::string thickSeparator = "=============================================================";
  std::string thinSeparator  = "-------------------------------------------------------------";

  ParameterList paramList;
  paramList.set("verbosity", "none");
  if (xmlFileName != "")
    Teuchos::updateParametersFromXmlFileAndBroadcast(xmlFileName, Teuchos::Ptr<ParameterList>(&paramList), *comm);

  out << "Parameter list:" << std::endl << paramList << std::endl;

  // =========================================================================
  // The LOOP
  // =========================================================================
  std::vector<std::string> reuseTypes, reuseNames;
  reuseTypes.push_back("none"); reuseNames.push_back("none");
  reuseTypes.push_back("S");    reuseNames.push_back("smoothers");
  reuseTypes.push_back("tP");   reuseNames.push_back("tentative P");
  reuseTypes.push_back("RP");   reuseNames.push_back("smoothed P and R");
  reuseTypes.push_back("RAP");  reuseNames.push_back("coarse grids");

  const size_t numSteps = 3;

  high_resolution_clock::time_point tc;
  std::vector<duration<double>> setup_time(reuseTypes.size()*numSteps);
  std::vector<duration<double>> solve_time(reuseTypes.size()*numSteps);
  std::vector<int>              num_its   (reuseTypes.size()*numSteps);

  for (size_t k = 0; k < reuseTypes.size(); k++) {
    out << thickSeparator << " " << reuseTypes[k] << " " << thickSeparator << std::endl;

    paramList.set("reuse: type", reuseTypes[k]);

    RCP<Matrix>       A;
    RCP<const Map>    map;
    RCP<MultiVector>  coordinates, nullspace;

    ConstructData(matrixType, galeriList, lib, comm, A, map, coordinates, nullspace);
    A->SetMaxEigenvalueEstimate(-one);

    RCP<Vector> X = VectorFactory::Build(map);
    X->setSeed(846930886);
    X->randomize();
    RCP<Vector> B = VectorFactory::Build(map);
    A->apply(*X, *B);

    RCP<Hierarchy> H = MueLu::CreateXpetraPreconditioner(A, paramList, coordinates);

    for (size_t t = 1; t < numSteps; t++) {
      out << thinSeparator << " Step " << t << " " << thinSeparator << std::endl;
      galeriList.set("stretchx", 1.0 + 10*t);
      galeriList.set("stretchy", 1.0);
      galeriList.set("stretchz", 0.5 / t);

      ConstructData(matrixType, galeriList, lib, comm, A, map, coordinates, nullspace);
      A->SetMaxEigenvalueEstimate(-one);

      tc = high_resolution_clock::now();
      if (solveType == "none")
        H = MueLu::CreateXpetraPreconditioner(A, paramList, coordinates);
      else
        MueLu::ReuseXpetraPreconditioner(A, H);
      setup_time[k*numSteps + t] = duration_cast<duration<double>>(high_resolution_clock::now() - tc);

      X->putScalar(zero);

      if (solveType == "none") {
        // Do nothing

      } else if (solveType == "standalone") {
        H->IsPreconditioner(false);

        tc = high_resolution_clock::now();
        H->Iterate(*B, *X, tol);
        solve_time[k*numSteps + t] = duration_cast<duration<double>>(high_resolution_clock::now() - tc);

      } else if (solveType == "cg" || solveType == "gmres") {
        H->IsPreconditioner(true);

#ifdef HAVE_MUELU_BELOS
        // Operator and Multivector type that will be used with Belos
        typedef MultiVector          MV;
        typedef Belos::OperatorT<MV> OP;

        // Define Operator and Preconditioner
        Teuchos::RCP<OP> belosOp   = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO>(A)); // Turns a Xpetra::Matrix object into a Belos operator
        Teuchos::RCP<OP> belosPrec = Teuchos::rcp(new Belos::MueLuOp <SC, LO, GO, NO>(H)); // Turns a MueLu::Hierarchy object into a Belos operator

        // Construct a Belos LinearProblem object
        RCP<Belos::LinearProblem<SC, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<SC, MV, OP>(belosOp, X, B));
        belosProblem->setRightPrec(belosPrec);

        bool set = belosProblem->setProblem();
        if (set == false) {
          out << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
          return EXIT_FAILURE;
        }

        // Belos parameter list
        Teuchos::ParameterList belosList;
        belosList.set("Maximum Iterations",    200); // Maximum number of iterations allowed
        belosList.set("Convergence Tolerance", tol);    // Relative convergence tolerance requested
        belosList.set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
        belosList.set("Output Frequency",      1);
        belosList.set("Output Style",          Belos::Brief);

        // Create an iterative solver manager
        RCP< Belos::SolverManager<SC, MV, OP> > solver;
        if (solveType == "cg") {
          solver = rcp(new Belos::PseudoBlockCGSolMgr   <SC, MV, OP>(belosProblem, rcp(&belosList, false)));
        } else if (solveType == "gmres") {
          solver = rcp(new Belos::BlockGmresSolMgr<SC, MV, OP>(belosProblem, rcp(&belosList, false)));
        }

        // Perform solve
        Belos::ReturnType ret = Belos::Unconverged;

        tc = high_resolution_clock::now();
        ret = solver->solve();
        solve_time[k*numSteps + t] = duration_cast<duration<double>>(high_resolution_clock::now() - tc);

        // Get the number of iterations for this solve.
        out << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;
        // Check convergence
        if (ret != Belos::Converged) {
          out << std::endl << "ERROR:  Belos did not converge! " << std::endl;
          num_its[k*numSteps+t] = -1;

        } else {
          out << std::endl << "SUCCESS:  Belos converged!" << std::endl;
          num_its[k*numSteps+t] = solver->getNumIters();
        }
#endif //ifdef HAVE_MUELU_BELOS
      } else {
        throw MueLu::Exceptions::RuntimeError("Unknown solver type: \"" + solveType + "\"");
      }

      out << "residual(A) = " << Utilities::ResidualNorm(*A, *X, *B)[0] << " [reuse \"" << reuseNames[k] << "\"]" << std::endl;
    }
  }
  for (size_t t = 1; t < numSteps; t++)
    for (size_t k = 0; k < reuseTypes.size(); k++)
      printf("step #%zu reuse \"%20s\": setup = %5.2le, solve = %5.2le, num_its = %3d\n", t, reuseNames[k].c_str(),
             setup_time[k*numSteps+t].count(), solve_time[k*numSteps+t].count(), num_its[k*numSteps+t]);

  return EXIT_SUCCESS;
}

int main(int argc, char* argv[]) {
  bool verbose = true;
  bool success = false;

  try {
    const bool throwExceptions = false;

    Teuchos::CommandLineProcessor clp(throwExceptions);
    Xpetra::Parameters xpetraParameters(clp);

    clp.recogniseAllOptions(false);
    switch (clp.parse(argc, argv, NULL)) {
      case Teuchos::CommandLineProcessor::PARSE_ERROR:               return EXIT_FAILURE;
      case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:
      case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION:
      case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
    }

    Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

    if (lib == Xpetra::UseEpetra) {
#ifdef HAVE_MUELU_EPETRA
      return main_<double,int,int,Xpetra::EpetraNode>(clp, argc, argv);
#else
      throw MueLu::Exceptions::RuntimeError("Epetra is not available");
#endif
    }

    if (lib == Xpetra::UseTpetra) {
      typedef KokkosClassic::DefaultNode::DefaultNodeType Node;

#ifndef HAVE_MUELU_EXPLICIT_INSTANTIATION
      return main_<double,int,long,Node>(clp, argc, argv);
#else
#  if defined(HAVE_MUELU_INST_DOUBLE_INT_INT)
      return main_<double,int,int,Node> (clp, argc, argv);
#elif defined(HAVE_MUELU_INST_DOUBLE_INT_LONGINT)
      return main_<double,int,long,Node>(clp, argc, argv);
#elif defined(HAVE_MUELU_INST_DOUBLE_INT_LONGLONGINT)
      return main_<double,int,long long,Node>(clp, argc, argv);
#else
      throw std::runtime_error("Found no suitable instantiation");
#endif
#endif
    }

  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
