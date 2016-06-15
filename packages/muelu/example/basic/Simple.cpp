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

#include <iostream>

// Tpetra provides distributed sparse linear algebra
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

// Belos provides Krylov solvers
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosPseudoBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosTpetraAdapter.hpp>

// MueLu main header: include most common header files in one line
#include <MueLu.hpp>
#include <MueLu_TpetraOperator.hpp>
#include <MueLu_CreateTpetraPreconditioner.hpp>

int main(int argc, char *argv[]) {

  // Define default types
  typedef double                                      scalar_type;
  typedef int                                         local_ordinal_type;
  typedef int                                         global_ordinal_type;
  typedef KokkosClassic::DefaultNode::DefaultNodeType node_type;

  // Convenient typedef's
  typedef Tpetra::Operator<scalar_type,local_ordinal_type,global_ordinal_type,node_type>    operator_type;
  typedef Tpetra::CrsMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type>   crs_matrix_type;
  typedef Tpetra::Vector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>      vector_type;
  typedef Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type> multivector_type;
  typedef Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type>                     driver_map_type;

  typedef MueLu::TpetraOperator<scalar_type, local_ordinal_type, global_ordinal_type, node_type> muelu_tpetra_operator_type;

  typedef Belos::LinearProblem<scalar_type, multivector_type, operator_type> linear_problem_type;
  typedef Belos::SolverManager<scalar_type, multivector_type, operator_type> belos_solver_manager_type;
  typedef Belos::PseudoBlockCGSolMgr<scalar_type, multivector_type, operator_type> belos_pseudocg_manager_type;
  typedef Belos::BlockGmresSolMgr<scalar_type, multivector_type, operator_type> belos_gmres_manager_type;

  //MueLu_UseShortNames.hpp wants these typedefs.
  typedef scalar_type         Scalar;
  typedef local_ordinal_type  LocalOrdinal;
  typedef global_ordinal_type GlobalOrdinal;
  typedef node_type           Node;
# include <MueLu_UseShortNames.hpp>

  using Teuchos::RCP; // reference count pointers
  using Teuchos::rcp; // reference count pointers
  Tpetra::global_size_t INVALID_GO = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();

  //
  // MPI initialization using Teuchos
  //

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  Teuchos::CommandLineProcessor clp(false);

  // Parameters

  global_ordinal_type numGlobalElements = 10000;
  global_ordinal_type maxIts            = 25;
  scalar_type tol                       = 1e-8;
  std::string solverOptionsFile         = "amg.xml";
  std::string krylovSolverType          = "cg";

  clp.setOption("nx",         &numGlobalElements, "number of matrix rows");
  clp.setOption("xmlFile",    &solverOptionsFile, "XML file containing MueLu solver parameters");
  clp.setOption("maxits",     &maxIts,            "maximum number of Krylov iterations");
  clp.setOption("tol",        &tol,               "tolerance for Krylov solver");
  clp.setOption("krylovType", &krylovSolverType,  "cg or gmres solver");

  switch (clp.parse(argc, argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:          break;
  }
  ParameterList mueluParams;
  Teuchos::updateParametersFromXmlFile(solverOptionsFile, Teuchos::inoutArg(mueluParams));

  //
  // Construct the problem
  //

  // Construct a Map that puts approximately the same number of equations on each processor
  RCP<const driver_map_type> map = rcp(new driver_map_type(INVALID_GO, numGlobalElements, 0, comm));
  

  // Get update list and number of local equations from newly created map.
  const size_t numMyElements = map->getNodeNumElements();
  Teuchos::ArrayView<const global_ordinal_type> myGlobalElements = map->getNodeElementList();

  // Create a CrsMatrix using the map, with a dynamic allocation of 3 entries per row
  RCP<crs_matrix_type> A = rcp(new crs_matrix_type(map, 3));

  // Add rows one-at-a-time
  for (size_t i = 0; i < numMyElements; i++) {
    if (myGlobalElements[i] == 0) {
      A->insertGlobalValues(myGlobalElements[i],
                            Teuchos::tuple<global_ordinal_type>(myGlobalElements[i], myGlobalElements[i] +1),
                            Teuchos::tuple<scalar_type> (2.0, -1.0));
    }
    else if (myGlobalElements[i] == numGlobalElements - 1) {
      A->insertGlobalValues(myGlobalElements[i],
                            Teuchos::tuple<global_ordinal_type>(myGlobalElements[i] -1, myGlobalElements[i]),
                            Teuchos::tuple<scalar_type> (-1.0, 2.0));
    }
    else {
      A->insertGlobalValues(myGlobalElements[i],
                            Teuchos::tuple<global_ordinal_type>(myGlobalElements[i] -1, myGlobalElements[i], myGlobalElements[i] +1),
                            Teuchos::tuple<scalar_type> (-1.0, 2.0, -1.0));
    }
  }

  // Complete the fill, ask that storage be reallocated and optimized
  A->fillComplete();

  //
  // Construct a multigrid preconditioner
  //

  // Multigrid Hierarchy
  RCP<muelu_tpetra_operator_type> M = MueLu::CreateTpetraPreconditioner((RCP<operator_type>)A, mueluParams);

  //
  // Set up linear problem Ax = b and associate preconditioner with it.
  //

  RCP<multivector_type> X = rcp(new multivector_type(map,1));
  RCP<multivector_type> B = rcp(new multivector_type(map,1));

  X->putScalar((scalar_type) 0.0);
  B->randomize();

  RCP<linear_problem_type> Problem = rcp(new linear_problem_type(A, X, B));
  Problem->setRightPrec(M);
  Problem->setProblem();

  //
  // Set up Krylov solver and iterate.
  //
  RCP<ParameterList> belosList = rcp(new ParameterList());
  belosList->set("Maximum Iterations",    maxIts); // Maximum number of iterations allowed
  belosList->set("Convergence Tolerance", tol);    // Relative convergence tolerance requested
  belosList->set("Verbosity",             Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
  belosList->set("Output Frequency",      1);
  belosList->set("Output Style",          Belos::Brief);
  belosList->set("Implicit Residual Scaling", "None");
  RCP<belos_solver_manager_type> solver;
  if (krylovSolverType == "cg")
    solver = rcp(new belos_pseudocg_manager_type(Problem, belosList));
  else if (krylovSolverType == "gmres")
    solver = rcp(new belos_gmres_manager_type(Problem, belosList));
  else
    throw std::invalid_argument("bad Krylov solver type");

  solver->solve();
  int numIterations = solver->getNumIters();

  Teuchos::Array<typename Teuchos::ScalarTraits<scalar_type>::magnitudeType> normVec(1);
  multivector_type Ax(B->getMap(),1);
  multivector_type residual(B->getMap(),1);
  A->apply(*X, residual);
  residual.update(1.0, *B, -1.0);
  residual.norm2(normVec);
  if (comm->getRank() == 0) {
    std::cout << "number of iterations = " << numIterations << std::endl;
    std::cout << "||Residual|| = " << normVec[0] << std::endl;
  }

  return EXIT_SUCCESS;
}
