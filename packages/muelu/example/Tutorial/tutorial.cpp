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

/*!
\example tutorial.cpp
\brief Driver examples for an application that wants to use MueLu as a preconditioner.
*/

//
// This example shows how to create a problem, set up a MueLu hierarchy with
// the help of an XML file, and to solve it using MueLu as a preconditioner.
//

#include <iostream>

// Define default data types
#include <Kokkos_DefaultNode.hpp>
#include <Kokkos_DefaultKernels.hpp>
typedef double                                                              Scalar;
typedef int                                                                 LocalOrdinal;
typedef int                                                                 GlobalOrdinal;
typedef KokkosClassic::DefaultNode::DefaultNodeType                         Node;
typedef KokkosClassic::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps  LocalMatOps;

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>

#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_Utilities.hpp>

#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>
#include <BelosMueLuAdapter.hpp>

int main(int argc, char *argv[]) {
  // Most MueLu and Xpetra classes are templated on some or all of the
  // following template types: Scalar, LocalOrdinal, GlobalOrdinal, Node, and
  // LocalMatOps. In order to make more concise, MueLu provides typedefs for
  // its classes, called short names. Thus, instead of writing, for instance,
  //
  //   MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> H;
  //
  // one can simply write
  //
  //   Hierarchy H;
  //
  // provided the MueLu_UseShortNames.hpp header and the appropriate class are
  // included.
#include <MueLu_UseShortNames.hpp>

  // These "using" declarations make the code more concise, in that you don't
  // have to write the namespace along with the class or object name. This is
  // especially helpful with commonly used things like std::endl or
  // Teuchos::RCP.
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::ParameterList;

  // Start up MPI, if using MPI. Trilinos doesn't have to be built
  // with MPI; it's called a "serial" build if you build without MPI.
  // GlobalMPISession hides this implementation detail.
  //
  // Note the third argument. If you pass GlobalMPISession the
  // address of an std::ostream, it will print a one-line status
  // message with the rank on each MPI process. This may be
  // undesirable if running with a large number of MPI processes.
  // You can avoid printing anything here by passing in either
  // NULL or the address of a Teuchos::oblackholestream.
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);

  // Get a pointer to the communicator object representing
  // MPI_COMM_WORLD.
  //
  // If we didn't build with MPI, we'll get a serial "communicator", i.e. a
  // communicator with size 1, whose only process has rank 0.
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // Get my process' rank. Equivalent to MPI_Comm_rank.
  const int myRank = comm->getRank();

  // Teuchos::ScalarTraits provides a traits interface to general scalar types.
  // This allows the same code to be valid for different scalar types, like
  // double, complex, or high precision (i.e., double-double)
  typedef Teuchos::ScalarTraits<SC> STS;
  SC zero = STS::zero(), one = STS::one();

  // Make an output stream (for verbose output) that only prints on
  // Proc 0 of the communicator.
  Teuchos::oblackholestream blackHole;
  std::ostream& out = (myRank == 0) ? std::cout : blackHole;

  // Teuchos provides an interface to get arguments from command line. The first
  // argument indicates that we don't want any exceptions to be thrown during
  // command line parsing
  Teuchos::CommandLineProcessor clp(false);

  // The list of valid parameters for the driver comes from three sources:
  //  - Galeri parameters
  //    These are responsible for constructing the problem data (i.e., matrix and coordinates).
  //  - Xpetra parameters
  //    Allow switching between Epetra and Tpetra
  //  - Driver parameters
  //    These are specified below
  GO nx = 100, ny = 100, nz = 100;
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace1D"); // manage parameters of the test case
  Xpetra::Parameters             xpetraParameters(clp);                          // manage parameters of Xpetra

  std::string xmlFileName = "tutorial1a.xml";   clp.setOption("xml",                   &xmlFileName,     "read parameters from a file [default = 'scalingTest.xml']");

  std::string solveType    = "cg";              clp.setOption("solver",                &solveType,       "solve type: (cg | gmres)");
  double      tol          = 1e-12;             clp.setOption("tol",                   &tol,             "solver convergence tolerance");
  int         maxIts       = 200;               clp.setOption("its",                   &maxIts,           "maximum number of solver iterations");

  std::string mapFile;                          clp.setOption("map",                   &mapFile,          "map data file");
  std::string matrixFile;                       clp.setOption("matrix",                &matrixFile,       "matrix data file");
  std::string coordFile;                        clp.setOption("coords",                &coordFile,        "coordinates data file");

  // Command line processor parsing stage
  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
  }

  // Retrieve matrix parameters (they may have been changed on the command line)
  ParameterList galeriList = galeriParameters.GetParameterList();

  // Construct the problem data
  // Typically, we construct matrix, coordinates, and nullspace, though only matrix is mandatory
  RCP<Matrix>      A;
  RCP<const Map>   map;
  RCP<MultiVector> coordinates;
  if (matrixFile.empty()) {
    // Use Galeri to construct the data
    // Galeri provides several pre-defined problem types, including some
    // stencil based ones (like Laplace matrices in 2D and 3D), and some others
    // (like elasticity in 2D and 3D)
    out << "========================================================\n" << xpetraParameters << galeriParameters;

    // Galeri will attempt to create a square-as-possible distribution of
    // subdomains di, e.g.,
    //                                 d1  d2  d3
    //                                 d4  d5  d6
    //                                 d7  d8  d9
    //                                 d10 d11 d12
    // A perfect distribution is only possible when the #processors is a
    // perfect square.  This *will* result in "strip" distribution if the
    // #processors is a prime number or if the factors are very different in
    // size. For example, np=14 will give a 7-by-2 distribution.  If you don't
    // want Galeri to do this, specify mx or my on the galeriList.
    std::string matrixType = galeriParameters.GetMatrixType();

    // Create map (Map) and coordinates (MultiVector)
    // Xpetra's Map and MultiVector follow Tpetra's interface
    if (matrixType == "Laplace1D") {
      map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian1D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("1D", map, galeriList);

    } else if (matrixType == "Laplace2D") {
      map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian2D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("2D", map, galeriList);

    } else if (matrixType == "Laplace3D") {
      map = Galeri::Xpetra::CreateMap<LO, GO, Node>(xpetraParameters.GetLib(), "Cartesian3D", comm, galeriList);
      coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<SC,LO,GO,Map,MultiVector>("3D", map, galeriList);
    }

    out << "Processor subdomains in x direction: " << galeriList.get<int>("mx") << std::endl
        << "Processor subdomains in y direction: " << galeriList.get<int>("my") << std::endl
        << "Processor subdomains in z direction: " << galeriList.get<int>("mz") << std::endl
        << "========================================================" << std::endl;

    // Construct the matrix based on the problem name and provided parameters
    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
        Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(galeriParameters.GetMatrixType(), map, galeriList);
    A = Pr->BuildMatrix();

  } else {
    // Read in data
    // This is a convenient way to use the driver, provided you have some
    // custom data you want to run.  Typically, you need a map (though you may
    // avoid that in a serial run), a matrix (in a MatrixMarket format), and
    // a file with coordinates
    RCP<const Map> map = (mapFile.empty() ? Teuchos::null : Utils2::ReadMap(mapFile, xpetraParameters.GetLib(), comm));

    A = Utils::Read(matrixFile, map);

    coordinates = Utils2::ReadMultiVector(coordFile, map);
  }

  // For scalar equations, we assume that the constant vector is a good
  // approximate nullspace Note, though, that we may skip this step as this is
  // also a default nullspace in MueLu, if user has not provided it.
  RCP<MultiVector> nullspace = MultiVectorFactory::Build(map, 1);
  nullspace->putScalar(one);

  out << "Galeri complete.\n========================================================" << std::endl;

  // Hierarchy is a key concept in MueLu. Once set up, it has all the data
  // required to solve Bx=y, where B is the preconditioner. MueLu provides
  // multiple interfaces to construct a Hierarchy. The first one is a C++
  // interface, which is suitable for advanced users. The second one is a XML
  // file interface, used in this driver. The documentation for the interface
  // and a list of parameter names is located in
  // trilinos/packages/muelu/doc/UsersGuide.
  //
  // HierarchyManager creates a Hierarchy object based on a provided XML file
  // parameters. We additionally provide the communicator object so that the
  // parameter list data can be redistributed to all processors.
  RCP<HierarchyManager> mueLuFactory = rcp(new ParameterListInterpreter(xmlFileName, *comm));

  // The CreateHierarchy() call creates a Hierarchy object. This object is not
  // set up, though, so it cannot be used yet.
  RCP<Hierarchy> H = mueLuFactory->CreateHierarchy();

  // All hierarchies have at least one level: the finest level, which must
  // contain users matrix. We use incrementing numbers for coarser levels. This
  // way, if we have L levels, level 0 corresponds to the finest level, and
  // level L-1 is the coarsest level.
  RCP<Level>     L = H->GetLevel(0);

  // Setting finest level matrix is mandatory in MueLu. However, setting
  // nullspace and coordinates is optional. If nullspace is not provided, it is
  // automatically generated to be a constant vector. Coordinates are needed
  // only in two cases: if one uses distance Laplacian filtering, or one uses
  // repartitioning (which is based on geometric partitioning).
  L->Set("A",           A);
  L->Set("Nullspace",   nullspace);
  L->Set("Coordinates", coordinates);

  // Now that we have set the mandatory data, we can construct the
  // preconditioner.  The result of this is a fully constructed hierarchy,
  // which can be used to solve the system.
  mueLuFactory->SetupHierarchy(*H);

  // Multigrid methods can be used in a standalon fashion, or they can be used to precondition a Krylov method.
  // Typically, the second method is preferred, which is indicated by the argument of IsPreconditioner function.
  H->IsPreconditioner(true);

  // To complete the problem description, we construct an initial guess X, and
  // a right hand size B. Both are base on the same map as is used for row map in A
  RCP<Vector> X = VectorFactory::Build(map);
  RCP<Vector> B = VectorFactory::Build(map);

  {
    // We set seed for reproducibility
    Utils::SetRandomSeed(*comm);
    X->randomize();
    A->apply(*X, *B, Teuchos::NO_TRANS, one, zero);

    Teuchos::Array<STS::magnitudeType> norms(1);
    B->norm2(norms);
    B->scale(one/norms[0]);
    X->putScalar(zero);
  }

  typedef MultiVector          MV;
  typedef Belos::OperatorT<MV> OP;

  // Internally, Belos uses its own types. We use special wrappers to wrap
  // MueLu and Xpetra objects into Belos operators
  RCP<OP> belosOp   = rcp(new Belos::XpetraOp<SC, LO, GO, NO, LMO>(A));
  RCP<OP> belosPrec = rcp(new Belos::MueLuOp <SC, LO, GO, NO, LMO>(H));

  // Construct a Belos LinearProblem object. This is a complete problem
  // formulation. All the data necessary to solve the system is now inside
  // that problem.
  RCP< Belos::LinearProblem<SC, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<SC, MV, OP>(belosOp, X, B));
  belosProblem->setRightPrec(belosPrec);

  // Prepare the linear problem to solve the linear system that was already passed in.
  bool set = belosProblem->setProblem();
  if (set == false) {
    out << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
    return EXIT_FAILURE;
  }

  // As many packages in Trilinos, Belos provides a parameter list interface
  // This allows for an easy tweaking of solver behaviour, such as changing the
  // number of iterations, or verbosity. For a full list of Belos parameters,
  // please consult the Belos package documentation.
  ParameterList belosList;
  belosList.set("Maximum Iterations",    maxIts); // Maximum number of iterations allowed
  belosList.set("Convergence Tolerance", tol);    // Relative convergence tolerance requested
  belosList.set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
  belosList.set("Output Frequency",      1);
  belosList.set("Output Style",          Belos::Brief);

  // Create an iterative solver manager
  // Belos provides multiple Krylov methods, with CG and GMRES being most
  // popular. Based on a command line argument, we select the proper one.  Note
  // that the C++ objects are templated on scalar, multivector and operator.
  // This is due to the traits interface in Belos.
  RCP<Belos::SolverManager<SC, MV, OP> > solver;
  if (solveType == "cg")
    solver = rcp(new Belos::PseudoBlockCGSolMgr<SC,MV,OP>(belosProblem, rcp(&belosList, false)));
  else if (solveType == "gmres")
    solver = rcp(new Belos::BlockGmresSolMgr   <SC,MV,OP>(belosProblem, rcp(&belosList, false)));

  // Finally, perform the solve. We wrap it in try-catch as Belos indicates
  // erros by throwing exceptions.
  Belos::ReturnType ret = Belos::Unconverged;
  try {
    ret = solver->solve();

    // Get the number of iterations for this solve.
    out << "Converged in " << solver->getNumIters() << " iterations" << std::endl;

  } catch(...) {
    out << std::endl << "ERROR:  Belos threw an error! " << std::endl;
  }

  // Check convergence
  if (ret != Belos::Converged)
    out << std::endl << "ERROR:  Belos did not converge! " << std::endl;
  else
    out << std::endl << "SUCCESS:  Belos converged!" << std::endl;

  // GlobalMPISession calls MPI_Finalize() in its destructor, if
  // appropriate. You don't have to do anything here!  Just return
  // from main(). Isn't that helpful?
  return 0;
}
