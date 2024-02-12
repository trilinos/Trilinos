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
//
// To compile and run this example, Trilinos must be configured with
// Tpetra, Amesos2, MueLu, Ifpack2, and Belos.
//
// This example will only work with Tpetra, not Epetra.
//
// Commonly used options are
//   --matrixType
//   --nx
//   --ny
//   --xmlFile
//
// "./MueLu_Repartition_ADR.exe --help" shows all supported command line options.
//

#include <iostream>

// Tpetra provides distributed sparse linear algebra
//#include <Tpetra_CrsMatrix.hpp>
//#include <Tpetra_Vector.hpp>

// Belos provides Krylov solvers
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosBiCGStabSolMgr.hpp>
#include <BelosTpetraAdapter.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>

// ADR subdirectory
#include "CreateADRMatrix.hpp"
#include <CreateADRMatrix.hpp>
#include <ADRProblemFactory.hpp>
#include <ADR_XpetraParameters.hpp>

// MueLu main header: include most common header files in one line
#include <MueLu.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>
#include <MueLu_Utilities.hpp>

#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_ExportFactory.hpp>
#include <Xpetra_MatrixFactory.hpp>

int main(int argc, char *argv[]) {
  // Define default types
  typedef double scalar_type;
  typedef int local_ordinal_type;
  typedef int global_ordinal_type;
  typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType node_type;

  // Convenient typedef's
  typedef Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type> operator_type;
  typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type, node_type> crs_matrix_type;
  typedef Tpetra::Vector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> vector_type;
  typedef Tpetra::MultiVector<scalar_type, local_ordinal_type, global_ordinal_type, node_type> multivector_type;
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> driver_map_type;

  typedef MueLu::TpetraOperator<scalar_type, local_ordinal_type, global_ordinal_type, node_type> muelu_tpetra_operator_type;
  typedef MueLu::Utilities<scalar_type, local_ordinal_type, global_ordinal_type, node_type> MueLuUtilities;

  typedef Belos::LinearProblem<scalar_type, multivector_type, operator_type> linear_problem_type;
  typedef Belos::SolverManager<scalar_type, multivector_type, operator_type> belos_solver_manager_type;
  typedef Belos::PseudoBlockCGSolMgr<scalar_type, multivector_type, operator_type> belos_pseudocg_manager_type;
  typedef Belos::BlockGmresSolMgr<scalar_type, multivector_type, operator_type> belos_gmres_manager_type;
  typedef Belos::BiCGStabSolMgr<scalar_type, multivector_type, operator_type> belos_bicgstab_manager_type;

  // MueLu_UseShortNames.hpp wants these typedefs.
  typedef scalar_type Scalar;
  typedef local_ordinal_type LocalOrdinal;
  typedef global_ordinal_type GlobalOrdinal;
  typedef node_type Node;
#include <MueLu_UseShortNames.hpp>

  typedef Galeri::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> GaleriXpetraProblem;
  typedef ADR::Xpetra::Problem<Map, CrsMatrixWrap, MultiVector> ADRXpetraProblem;

  using Teuchos::RCP;  // reference count pointers
  using Teuchos::rcp;  // reference count pointers

  //
  // MPI initialization using Teuchos
  //

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int mypid                           = comm->getRank();
  /*
    int subCommRank[3]={0,1,2};
    Teuchos::ArrayView<int> arraySubCommRank(subCommRank, 3);
    auto subComm = comm->createSubcommunicator(arraySubCommRank);
  */
  Teuchos::CommandLineProcessor clp(false);

  // Problem definition
  std::string problem_type = "ADR2D";

  // Parameters

  Scalar Lx = 6.0;
  Scalar Ly = 6.0;
  Scalar Lz = 6.0;

  global_ordinal_type nx = 50;
  global_ordinal_type ny = 50;
  global_ordinal_type nz = 50;

  int keep_boundary = 0;
  Scalar stretchx   = (Scalar)Lx / nx;
  Scalar stretchy   = (Scalar)Ly / ny;
  Scalar stretchz   = (Scalar)Lz / nz;

  ADR::Xpetra::Parameters<GO> matrixParameters(clp, nx, ny, nz, problem_type, keep_boundary, stretchx, stretchy, stretchz);  // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);                                                                                  // manage parameters of xpetra

  global_ordinal_type maxIts    = 10000;
  scalar_type tol               = 1e-10;
  std::string solverOptionsFile = "dd.xml";
  std::string krylovSolverType  = "bicgstab";

  clp.setOption("xmlFile", &solverOptionsFile, "XML file containing MueLu solver parameters");
  clp.setOption("maxits", &maxIts, "maximum number of Krylov iterations");
  clp.setOption("tol", &tol, "tolerance for Krylov solver");
  clp.setOption("krylovType", &krylovSolverType, "cg or gmres solver");

  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED: return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL: break;
  }

  if (xpetraParameters.GetLib() == Xpetra::UseEpetra) {
    throw std::invalid_argument("This example only supports Tpetra.");
  }

  ParameterList mueluParams;
  Teuchos::updateParametersFromXmlFile(solverOptionsFile, Teuchos::inoutArg(mueluParams));

  //
  // Construct the problem
  //

  global_ordinal_type indexBase = 0;
  RCP<const Map> xpetraMap;
  std::vector<global_ordinal_type> indices;

  // Creation of the map where processor 0 gets nothing at the fine level
  if (comm->getSize() > 1) {
    indices.reserve(static_cast<int>(matrixParameters.GetNumGlobalElements() / (comm->getSize() - 1) + 1));
    if (mypid != 0 && mypid != comm->getSize() - 1)
      for (int i = 0; i <= (static_cast<int>(matrixParameters.GetNumGlobalElements() / (comm->getSize() - 1))) - 1; ++i)
        indices.emplace_back((mypid - 1) * static_cast<int>(matrixParameters.GetNumGlobalElements() / (comm->getSize() - 1)) + i);

    if (mypid == comm->getSize() - 1)
      for (int i = (mypid - 1) * static_cast<int>(matrixParameters.GetNumGlobalElements() / (comm->getSize() - 1)); i != matrixParameters.GetNumGlobalElements(); ++i)
        indices.emplace_back(i);

    indices.shrink_to_fit();

    Teuchos::ArrayView<const global_ordinal_type> elementList(indices);
    xpetraMap = MapFactory::Build(Xpetra::UseTpetra, matrixParameters.GetNumGlobalElements(), elementList, indexBase, comm);
  } else if (comm->getSize() == 1) {
    xpetraMap = MapFactory::Build(Xpetra::UseTpetra, matrixParameters.GetNumGlobalElements(), indexBase, comm);
  }

  //============================================================================================================

  RCP<MultiVector> coordinates;

  // RCP<const Map> xpetraMap = MapFactory::Build(Xpetra::UseTpetra, matrixParameters.GetNumGlobalElements(), indexBase, comm);
  // RCP<const Map> xpetraMap = Galeri::Xpetra::CreateMap<local_ordinal_type, global_ordinal_type, node_type>(Xpetra::UseTpetra, "Cartesian2D", comm, matrixParameters.GetParameterList());

  std::cout << "Before coordinates" << std::endl;
  coordinates = Galeri::Xpetra::Utils::CreateCartesianCoordinates<scalar_type, local_ordinal_type, global_ordinal_type, Map, MultiVector>("2D", xpetraMap, matrixParameters.GetParameterList());
  std::cout << "After coordinates" << std::endl;

  RCP<ADRXpetraProblem> Pr = ADR::Xpetra::BuildProblem<scalar_type, local_ordinal_type, global_ordinal_type, Map, CrsMatrixWrap, MultiVector>(matrixParameters.GetMatrixType(), xpetraMap, matrixParameters.GetParameterList());
  RCP<Matrix> xpetraA      = Pr->BuildMatrix();

  RCP<crs_matrix_type> A         = MueLuUtilities::Op2NonConstTpetraCrs(xpetraA);
  RCP<const driver_map_type> map = MueLuUtilities::Map2TpetraMap(*xpetraMap);

  //
  // Construct a multigrid preconditioner
  //

  RCP<muelu_tpetra_operator_type> M = MueLu::CreateTpetraPreconditioner((RCP<operator_type>)A, mueluParams, Utilities::MV2NonConstTpetraMV(coordinates));

  RCP<multivector_type> X = rcp(new multivector_type(map, 1));
  RCP<multivector_type> B = rcp(new multivector_type(map, 1));

  X->putScalar((scalar_type)0.0);
  B->randomize();

  RCP<linear_problem_type> Problem = rcp(new linear_problem_type(A, X, B));
  Problem->setRightPrec(M);
  Problem->setProblem();

  //
  // Set up Krylov solver and iterate.
  //

  RCP<ParameterList> belosList = rcp(new ParameterList());
  belosList->set("Maximum Iterations", maxIts);  // Maximum number of iterations allowed
  belosList->set("Convergence Tolerance", tol);  // Relative convergence tolerance requested
  belosList->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
  belosList->set("Output Frequency", 1);
  belosList->set("Output Style", Belos::Brief);
  belosList->set("Implicit Residual Scaling", "None");
  RCP<belos_solver_manager_type> solver;
  if (krylovSolverType == "cg")
    solver = rcp(new belos_pseudocg_manager_type(Problem, belosList));
  else if (krylovSolverType == "gmres")
    solver = rcp(new belos_gmres_manager_type(Problem, belosList));
  else if (krylovSolverType == "bicgstab")
    solver = rcp(new belos_bicgstab_manager_type(Problem, belosList));
  else
    throw std::invalid_argument("bad Krylov solver type");

  solver->solve();
  int numIterations = solver->getNumIters();

  Teuchos::Array<typename Teuchos::ScalarTraits<scalar_type>::magnitudeType> normVec(1);
  multivector_type Ax(B->getMap(), 1);
  multivector_type residual(B->getMap(), 1);
  A->apply(*X, residual);
  residual.update(1.0, *B, -1.0);
  residual.norm2(normVec);
  if (mypid == 0) {
    std::cout << "number of iterations = " << numIterations << std::endl;
    std::cout << "||Residual|| = " << normVec[0] << std::endl;
  }

  return EXIT_SUCCESS;
}
