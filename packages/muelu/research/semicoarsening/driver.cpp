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

This example shows how to create a problem, set up a MueLu hierarchy
with the help of an XML file, and to solve it using MueLu as a
preconditioner.  It is meant to be used with one of the .xml files in
this directory.
*/

//
// This example shows how to create a problem, set up a MueLu hierarchy with
// the help of an XML file, and to solve it using MueLu as a preconditioner.
//

#include <iostream>

// Define default data types
//#include <Kokkos_DefaultNode.hpp>

//typedef double                                                              Scalar;
//typedef int                                                                 LocalOrdinal;
//typedef int                                                                 GlobalOrdinal;
//typedef KokkosClassic::DefaultNode::DefaultNodeType                         Node;

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_IO.hpp>

#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <MueLu_Utilities.hpp>
#include <MueLu_SemiCoarsenPFactory_decl.hpp>
#include "Teuchos_ScalarTraits.hpp"

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>
#include <BelosMueLuAdapter.hpp>
#endif

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
int main_(Teuchos::CommandLineProcessor &clp, Xpetra::UnderlyingLib lib, int argc, char *argv[]) {
  // Most MueLu and Xpetra classes are templated on some or all of the
  // following template types: Scalar, LocalOrdinal, GlobalOrdinal,
  // and Node. In order to make types more concise, MueLu has an
  // option to let users refer to its classes by typedefs, which we
  // call "short names."  Thus, instead of writing, for instance,
  //
  //   MueLu::Hierarchy<Scalar,LocalOrdinal,GlobalOrdinal,Node> H;
  //
  // one can simply write
  //
  //   Hierarchy H;
  //
  // provided the MueLu_UseShortNames.hpp header and the appropriate
  // class are included.  If one wishes to use short names, one must
  // include this header file before referring to any of the short
  // names.  One may include it anywhere in the program, but it is
  // better to include it in the scope where one wants the typedefs
  // defined.  It's better not to include it at global scope, in order
  // to avoid polluting the global scope with typedefs.
#include <MueLu_UseShortNames.hpp>

  // These "using" declarations make the code more concise, in that
  // you don't have to write the namespace along with the class or
  // object name. This is especially helpful with commonly used things
  // like std::endl or Teuchos::RCP.
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

  // Get a pointer to the communicator object representing
  // MPI_COMM_WORLD.
  //
  // If we didn't build with MPI, we'll get a serial "communicator",
  // i.e. a communicator with size 1, whose only process has rank 0.
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // Get my process' rank. Equivalent to MPI_Comm_rank.
  const int myRank = comm->getRank();

  // Teuchos::ScalarTraits provides a traits interface to general
  // scalar types.  This allows the same code to be valid for
  // different scalar types, like double, complex, or high precision
  // (i.e., double-double).
  typedef Teuchos::ScalarTraits<SC> STS;
  typedef typename STS::magnitudeType MT;
  SC zero = STS::zero(), one = STS::one();

  // Make an output stream (for verbose output) that only prints on
  // Process 0 of the communicator.
  Teuchos::oblackholestream blackHole;
  std::ostream& out = (myRank == 0) ? std::cout : blackHole;

  // Teuchos provides an interface to get arguments from the command
  // line. The first argument indicates that we don't want any
  // exceptions to be thrown during command line parsing.

  // The list of valid parameters for the driver comes from three sources:
  //
  //  - Galeri parameters: These are responsible for constructing the
  //    problem data (i.e., matrix and coordinates).
  //  - Xpetra parameters: Allow switching between Epetra and Tpetra.
  //  - Driver parameters: These are specified below.
  GO nx = 6, ny = 5, nz = 163;
  Galeri::Xpetra::Parameters<GO> galeriParameters(clp, nx, ny, nz, "Laplace3D"); // manage parameters of the test case
  Xpetra::Parameters             xpetraParameters(clp);                          // manage parameters of Xpetra

  std::string xmlFileName = "driver1.xml";   clp.setOption("xml",                   &xmlFileName,     "read parameters from a file [default = 'scalingTest.xml']");

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

  // Retrieve matrix parameters (they may have been changed on the
  // command line).
  ParameterList galeriList = galeriParameters.GetParameterList();
//galeriList.set("mx",6); galeriList.set("my",5); galeriList.set("mz",1);

  // Construct the problem data.  Typically, we construct the matrix,
  // coordinates, and nullspace, though only the matrix is mandatory.
  RCP<Matrix>      A;
  RCP<const Map>   map;
  RCP<MultiVector> coordinates;
  if (matrixFile.empty()) {
    // Use Galeri to construct the data
    //
    // Galeri provides several pre-defined problem types, including
    // some stencil based ones (like Laplace matrices in 2D and 3D),
    // and some others (like elasticity in 2D and 3D)
    out << "========================================================\n" << xpetraParameters << galeriParameters;

    // Galeri will attempt to create a square-as-possible distribution of
    // subdomains di, e.g.,
    //
    //                                 d1  d2  d3
    //                                 d4  d5  d6
    //                                 d7  d8  d9
    //                                 d10 d11 d12
    //
    // A perfect distribution is only possible when the number of
    // processes is a perfect square.  This *will* result in "strip"
    // distribution if the number of processes is a prime number or if
    // the factors are very different in size. For example, np=14 will
    // give a 7-by-2 distribution.  If you don't want Galeri to do
    // this, specify mx or my on the galeriList.
    std::string matrixType = galeriParameters.GetMatrixType();

    // Create map (Map) and coordinates (MultiVector).  Xpetra's Map
    // and MultiVector imitate Tpetra's interface.
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

    // Construct the matrix based on the problem name and provided
    // parameters.
    RCP<Galeri::Xpetra::Problem<Map,CrsMatrixWrap,MultiVector> > Pr =
      Galeri::Xpetra::BuildProblem<SC,LO,GO,Map,CrsMatrixWrap,MultiVector>(galeriParameters.GetMatrixType(), map, galeriList);
    A = Pr->BuildMatrix();

  } else {
    // Read in data
    //
    // This is a convenient way to use the driver, provided you have
    // some custom data you want to run.  Typically, you need a Map
    // (though you may avoid that in a serial run), a matrix (in a
    // MatrixMarket format), and a file with coordinates.
    if (!mapFile.empty())
      map = Xpetra::IO<SC,LO,GO,Node>::ReadMap(mapFile, xpetraParameters.GetLib(), comm);

    A = Xpetra::IO<SC,LO,GO,Node>::Read(matrixFile, map);

    if (!coordFile.empty())
      coordinates = Xpetra::IO<SC,LO,GO,Node>::ReadMultiVector(coordFile, map);
  }

  // For scalar equations, we assume that the constant vector is a
  // good approximate nullspace.  Note, though, that we may skip this
  // step as this is also a default nullspace in MueLu, if user has
  // not provided it.
  RCP<MultiVector> nullspace = MultiVectorFactory::Build(map, 1);
  nullspace->putScalar(one);

  out << "Galeri complete.\n========================================================" << std::endl;
std::cout << coordinates << std::endl;

  // Hierarchy is a key concept in MueLu. Once set up, it has all the
  // data required to apply the preconditioner. MueLu provides
  // multiple interfaces to construct a Hierarchy. The first one is a
  // C++ interface, which is suitable for advanced users. The second
  // one is a XML file interface, used in this driver. You may find
  // the documentation for the interface and a list of parameter names
  // in Trilinos/packages/muelu/doc/UsersGuide.
  //
  // HierarchyManager creates a Hierarchy object based on an XML file
  // containing the parameters.  Only Process 0 reads from this file.
  // We additionally provide a communicator, so that Process 0 can
  // broadcast the parameters to all other processes.
  RCP<HierarchyManager> mueLuFactory =
      rcp(new ParameterListInterpreter(xmlFileName, *comm));

  // The CreateHierarchy() call creates a Hierarchy object. This
  // object is not set up, though, so it cannot be used yet.
  RCP<Hierarchy> H = mueLuFactory->CreateHierarchy();

  // All hierarchies have at least one level: the finest level, which
  // must contain the user's matrix.  The finest level is level 0, and
  // subsequent coarser levels have increasing indices.  This way, if
  // we have L levels, level 0 corresponds to the finest level, and
  // level L-1 is the coarsest level.
  RCP<Level>     L = H->GetLevel(0);

  // Setting the finest level matrix is mandatory in MueLu. However,
  // setting the nullspace and coordinates is optional. If nullspace
  // is not provided, it is automatically generated to be a constant
  // vector. Coordinates are needed only in two cases: if one uses
  // distance Laplacian filtering, or one uses repartitioning (which
  // is based on geometric partitioning).
  L->Set("A",         A);
  L->Set("Nullspace", nullspace);

/*
  if (!coordinates.is_null())
    L->Set("Coordinates", coordinates);
*/

Teuchos::ArrayRCP<LO> SemiCoarsenInfo = Teuchos::arcp<LO>(3);
SemiCoarsenInfo[NUM_ZPTS] = nz;
SemiCoarsenInfo[ORIENTATION] = VERTICAL;
//L->Set("SemiCoarsenInfo",SemiCoarsenInfo, MueLu::NoFactory::get() );
L->Set("SemiCoarsenInfo",SemiCoarsenInfo);

printf("before level print\n");
L->print(std::cout,Teuchos::VERB_EXTREME);
printf("after level print\n");

  // Now that we have set the mandatory data, we can construct the
  // preconditioner.  The result of this is a fully constructed
  // hierarchy, which can be used to solve the linear system.
  mueLuFactory->SetupHierarchy(*H);

  // Multigrid methods can be used in a standalone fashion, or they
  // can be used to precondition a Krylov method.  Typically, the
  // second method is preferred, which is indicated by the argument of
  // IsPreconditioner function.
  H->IsPreconditioner(true);

  // To complete the problem description, we construct an initial
  // guess X and a right-hand side B (of the linear system AX=B). Both
  // have the same Map as A's row Map.
  RCP<Vector> X = VectorFactory::Build(map);
  RCP<Vector> B = VectorFactory::Build(map);

  // For this test problem, we determine the exact solution first (as
  // a random vector), and compute B so that AX=B exactly.  We then
  // make X an initial guess vector by filling it with zeros.  You may
  // use any initial guess you like.
  {
    // We set seed for reproducibility
    Utilities::SetRandomSeed(*comm);
    X->randomize();
    A->apply(*X, *B, Teuchos::NO_TRANS, one, zero);

    Teuchos::Array<MT> norms(1);
    B->norm2(norms);
    B->scale(one/norms[0]);
    X->putScalar(zero);
  }
#ifdef HAVE_MUELU_BELOS
  typedef MultiVector          MV;
  typedef Belos::OperatorT<MV> OP;

  // Internally, Belos uses its own types. We use special wrappers to
  // wrap MueLu and Xpetra objects into Belos operators.  This is only
  // one of many different ways one could choose to wrap MueLu and
  // Xpetra objects in order to get them to work with Belos.
  RCP<OP> belosOp   = rcp(new Belos::XpetraOp<SC, LO, GO, NO>(A));
  RCP<OP> belosPrec = rcp(new Belos::MueLuOp <SC, LO, GO, NO>(H));

  // Construct a Belos LinearProblem object. This is a complete
  // problem formulation. All the data necessary to solve the system
  // are now inside that problem.
  RCP< Belos::LinearProblem<SC, MV, OP> > belosProblem =
    rcp (new Belos::LinearProblem<SC, MV, OP> (belosOp, X, B));
  belosProblem->setRightPrec (belosPrec);

  // Prepare the linear problem to solve the linear system that was
  // already passed in.
  bool set = belosProblem->setProblem ();
  if (!set) {
    out << "\nERROR:  Belos::LinearProblem failed to set up correctly!" << std::endl;
    return EXIT_FAILURE;
  }

  // As many packages in Trilinos, Belos provides a parameter list
  // interface.  This lets users easily change solver behavior, like
  // the number of iterations or its verbosity. For a full list of
  // Belos parameters, please consult the Belos package documentation.
  ParameterList belosList;
  belosList.set("Maximum Iterations",    maxIts); // Maximum number of iterations allowed
  belosList.set("Convergence Tolerance", tol);    // Relative convergence tolerance requested
  belosList.set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
  belosList.set("Output Frequency",      1);
  belosList.set("Output Style",          Belos::Brief);

  // Create a Belos iterative linear solver
  //
  // Belos provides multiple iterative linear solvers, with CG and
  // GMRES being most popular. Based on a command line argument, we
  // select the proper one.  Note that the C++ objects are templated
  // on scalar, multivector and operator.  This is due to the traits
  // interface in Belos.
  RCP<Belos::SolverManager<SC, MV, OP> > solver;
  if (solveType == "cg")
    solver = rcp(new Belos::PseudoBlockCGSolMgr<SC,MV,OP>(belosProblem, rcp(&belosList, false)));
  else if (solveType == "gmres")
    solver = rcp(new Belos::BlockGmresSolMgr   <SC,MV,OP>(belosProblem, rcp(&belosList, false)));

  // Finally, perform the solve.  We wrap it in a try-catch block,
  // since Belos indicates erros by throwing exceptions.
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
#else
  out << std::endl << "MueLu has been compiled without Belos support!" << std::endl;
#endif
  // GlobalMPISession calls MPI_Finalize() in its destructor, if
  // appropriate. You don't have to do anything here!  Just return
  // from main(). Isn't that helpful?
  return 0;
}


//- -- --------------------------------------------------------
#define MUELU_AUTOMATIC_TEST_ETI_NAME main_
#include "MueLu_Test_ETI.hpp"

int main(int argc, char *argv[]) {
  return Automatic_Test_ETI(argc,argv);
}



