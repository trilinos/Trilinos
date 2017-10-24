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
#include <complex>

// Teuchos
#include <Teuchos_StandardCatchMacros.hpp>

// Tpetra
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <TpetraExt_MatrixMatrix.hpp>

// MueLu
#include <MueLu_RefMaxwell.hpp>
#include <MueLu_UseDefaultTypesComplex.hpp>
#include <MueLu_Exceptions.hpp>

// Belos
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>

int main(int argc, char *argv[]) {

#if defined(HAVE_MUELU_TPETRA) && defined(HAVE_MUELU_IFPACK2)

#include <MueLu_UseShortNames.hpp>

  typedef Tpetra::Map<LO,GO,NO>               TMap;
  typedef Tpetra::MultiVector<SC,LO,GO,NO>    TMV;
  typedef Tpetra::CrsMatrix<SC,LO,GO,NO>      TCRS;
  typedef Tpetra::Operator<SC,LO,GO,NO>       OP;
  typedef Belos::LinearProblem<SC,TMV,OP>     BelosProblem;
  typedef Belos::SolverManager<SC,TMV,OP>     BelosManager;
  typedef Belos::SolverFactory<SC,TMV,OP>     BelosFactory;

  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);

  bool success = false;
  bool verbose = true;
  try {
    RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    int commrank = comm->getRank();

    // Read matrices in from files
    Xpetra::global_size_t nedges=540, nnodes=216;
    // maps for nodal and edge matrices
    RCP<TMap> edge_map = rcp( new TMap(nedges,0,comm) );
    RCP<TMap> node_map = rcp( new TMap(nnodes,0,comm) );
    // edge stiffness matrix
    RCP<Matrix> S_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read("S.txt", edge_map);
    // edge mass matrix
    RCP<Matrix> M1_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read("M1.txt", edge_map);
    // nodal mass matrix
    RCP<Matrix> M0_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read("M0.txt", node_map);
    // gradient matrix
    RCP<Matrix> D0_Matrix = Xpetra::IO<SC, LO, GO, NO>::Read("D0.txt", edge_map, Teuchos::null, node_map, edge_map);
    // coordinates
    RCP<MultiVector> coords = Xpetra::IO<SC, LO, GO, NO>::ReadMultiVector("coords.txt", node_map);

    // build lumped mass matrix inverse (M0inv_Matrix)
    RCP<TMV> ones = rcp( new TMV(node_map,1) );
    RCP<TMV> diag = rcp( new TMV(node_map,1) );
    RCP<TMV> invdiag = rcp( new TMV(node_map,1) );
    ones->putScalar((SC)1.0);
    M0_Matrix->apply(*ones,*diag);
    invdiag->reciprocal(*diag);
    Teuchos::ArrayRCP<const SC> invdiags = invdiag->getData(0);
    GO row, col;
    SC entry;
    RCP<TCRS> M0inv_Matrix = rcp( new TCRS(node_map,1) );
    for(int i=0; i<nnodes; i++) {
      row = i;
      col = i;
      if(node_map->isNodeGlobalElement(i)) {
        LocalOrdinal lclidx = node_map->getLocalElement(i);
        std::complex<double> centry = invdiags[lclidx];
        M0inv_Matrix -> insertGlobalValues(row,
            Teuchos::ArrayView<LO>(&col,1),
            Teuchos::ArrayView<SC>(&centry,1));
      }
    }
    M0inv_Matrix->fillComplete();
    // build stiffness plus mass matrix (SM_Matrix)
    RCP<TCRS> SM_Matrix = rcp( new TCRS(edge_map,100) );
    std::complex<double> omega(0.0,2.0*M_PI);
    Tpetra::MatrixMatrix::Add(*S_Matrix,false,(SC)1.0,*M1_Matrix,false,omega,SM_Matrix);
    SM_Matrix->fillComplete();

    // set parameters
    Teuchos::ParameterList params, params11, params22;
    params.set("refmaxwell: disable add-on",false);
    params.set("refmaxwell: max coarse size",25);
    params.set("max levels",4);
    params11.set("smoother: type","KRYLOV");
    params11.set("smoother: type","KRYLOV");
    //    params11.set("krylov: number of iterations",3);
    //    params22.set("krylov: number of iterations",3);
    params.set("refmaxwell: 11list",params11);
    params.set("refmaxwell: 22list",params22);
    // construct preconditioner
    RCP<MueLu::RefMaxwell<SC,LO,GO,NO> > preconditioner
      = rcp( new MueLu::RefMaxwell<SC,LO,GO,NO>(SM_Matrix,D0_Matrix,M0inv_Matrix,
            M1_Matrix,Teuchos::null,coords,params) );

    // setup LHS, RHS
    RCP<TMV> vec = rcp( new TMV(edge_map,1) );
    vec -> putScalar((SC)1.0);
    RCP<TMV> B = rcp( new TMV(edge_map,1) );
    SM_Matrix->apply(*vec,*B);
    RCP<TMV> X = rcp( new TMV(edge_map,1) );
    X -> putScalar((SC)0.0);
    // Belos linear problem
    RCP<BelosProblem> problem = rcp( new BelosProblem() );
    problem -> setOperator( SM_Matrix );
    problem -> setRightPrec( preconditioner );
    problem -> setProblem( X, B );
    // Belos solver
    RCP<BelosManager> solver;
    RCP<BelosFactory> factory = rcp( new BelosFactory() );
    RCP<Teuchos::ParameterList> belosParams
      = rcp( new Teuchos::ParameterList() );
    belosParams->set("Maximum Iterations", 100);
    belosParams->set("Convergence Tolerance",1e-9);
    belosParams->set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
    belosParams->set("Output Frequency",1);
    belosParams->set("Output Style",Belos::Brief);
    solver = factory->create("Flexible GMRES",belosParams);
    // set problem and solve
    solver -> setProblem( problem );
    Belos::ReturnType status = solver -> solve();
    int iters = solver -> getNumIters();
    success = (iters<20 && status == Belos::Converged);
    if (commrank == 0) {
      if (success)
        std::cout << "SUCCESS! Belos converged in " << iters << " iterations." << std::endl;
      else
        std::cout << "FAILURE! Belos did not converge fast enough." << std::endl;
    }
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
#endif
} // main
