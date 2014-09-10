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
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int commrank = comm->getRank();

  // Read matrices in from files
  std::ifstream inputfile;
  int nnz_grad=1080, nnz_nodes=4096, nnz_edges=13440;
  int nedges=540, nnodes=216, row, col;
  double entry, x, y, z;
  // maps for nodal and edge matrices
  RCP<TMap> edge_map = rcp( new TMap(nedges,0,comm) );
  RCP<TMap> node_map = rcp( new TMap(nnodes,0,comm) );
  // edge stiffness matrix
  RCP<TCRS> S_Matrix = rcp( new TCRS(edge_map,100) );
  inputfile.open("S.txt");
  for(int i=0; i<nnz_edges; i++) {
    inputfile >> row >> col >> entry ;
    row=row-1;
    col=col-1;
    std::complex<double> centry(entry,0.0);
    if(edge_map->isNodeGlobalElement(row)) {
      S_Matrix->insertGlobalValues(row,
				   Teuchos::ArrayView<LO>(&col,1),
				   Teuchos::ArrayView<SC>(&centry,1));
    }
  }
  S_Matrix->fillComplete();
  inputfile.close();
  // edge mass matrix
  RCP<TCRS> M1_Matrix = rcp( new TCRS(edge_map,100) );
  inputfile.open("M1.txt");
  for(int i=0; i<nnz_edges; i++) {
    inputfile >> row >> col >> entry ;
    row=row-1;
    col=col-1;
    std::complex<double> centry(entry,0.0);
    if(edge_map->isNodeGlobalElement(row)) {
      M1_Matrix->insertGlobalValues(row,
				    Teuchos::ArrayView<LO>(&col,1),
				    Teuchos::ArrayView<SC>(&centry,1));
    }
  }
  M1_Matrix->fillComplete();
  inputfile.close();
  // nodal mass matrix
  RCP<TCRS> M0_Matrix = rcp( new TCRS(node_map,100) );
  inputfile.open("M0.txt");
  for(int i=0; i<nnz_nodes; i++) {
    inputfile >> row >> col >> entry ;
    row=row-1;
    col=col-1;
    std::complex<double> centry(entry,0.0);
    if(node_map->isNodeGlobalElement(row)) {
      M0_Matrix->insertGlobalValues(row,
				    Teuchos::ArrayView<LO>(&col,1),
				    Teuchos::ArrayView<SC>(&centry,1));
    }
  }
  M0_Matrix->fillComplete();
  inputfile.close();
  // gradient matrix
  RCP<TCRS> D0_Matrix = rcp( new TCRS(edge_map,2) );
  inputfile.open("D0.txt");
  for(int i=0; i<nnz_grad; i++) {
    inputfile >> row >> col >> entry ;
    row=row-1;
    col=col-1;
    std::complex<double> centry(entry,0.0);
    if(edge_map->isNodeGlobalElement(row)) {
      D0_Matrix->insertGlobalValues(row,
				    Teuchos::ArrayView<LO>(&col,1),
				    Teuchos::ArrayView<SC>(&centry,1));
    }
  }
  D0_Matrix->fillComplete(node_map,edge_map);
  inputfile.close();
  // coordinates
  RCP<TMV> coords = rcp( new TMV(node_map,3) );
  inputfile.open("coords.txt");
  for(int i=0; i<nnodes; i++) {
    inputfile >> x >> y >> z ;
    std::complex<double> cx(x,0.0), cy(y,0.0), cz(z,0.0);
    if(node_map->isNodeGlobalElement(i)) {
      coords->replaceGlobalValue(i,0,cx);
      coords->replaceGlobalValue(i,1,cy);
      coords->replaceGlobalValue(i,2,cz);
    }
  }
  inputfile.close();
  // build lumped mass matrix inverse (M0inv_Matrix)
  RCP<TMV> ones = rcp( new TMV(node_map,1) );
  RCP<TMV> diag = rcp( new TMV(node_map,1) );
  RCP<TMV> invdiag = rcp( new TMV(node_map,1) );
  ones->putScalar((SC)1.0);
  M0_Matrix->apply(*ones,*diag);
  invdiag->reciprocal(*diag);
  Teuchos::ArrayRCP<const SC> invdiags = invdiag->getData(0);
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
  params.set("refmaxwell: max levels",4);
  params.set("refmaxwell: edge smoother","KRYLOV");
  params.set("refmaxwell: node smoother","KRYLOV");
  params11.set("krylov: number of iterations",3);
  params22.set("krylov: number of iterations",3);
  params.set("refmaxwell: edge smoother list",params11);
  params.set("refmaxwell: node smoother list",params22);
  // construct preconditioner
  RCP<RefMaxwell> preconditioner
    = rcp( new RefMaxwell(SM_Matrix,D0_Matrix,M0inv_Matrix,
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
  if(iters<20 && status == Belos::Converged) {
    if(commrank==0) {
      std::cout<<"SUCCESS! Belos converged in "<<iters<<" iterations."<<std::endl;
    }
  }
  else {
    throw(MueLu::Exceptions::RuntimeError("FAIL! Belos did not converge fast enough."));
  }
  
#endif

} // main
