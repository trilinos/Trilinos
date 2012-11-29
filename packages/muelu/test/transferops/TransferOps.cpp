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
/*
 * TransferOpsTest.cpp
 *
 *  Created on: 08.10.2011
 *      Author: tobias
 */

#include <unistd.h>
#include <iostream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Parameters.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>

// MueLu
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Memory.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Exceptions.hpp"

// Belos stuff
#ifdef HAVE_MUELU_BELOS
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosMueLuAdapter.hpp" // this header defines Belos::MueLuOp()
#endif

//
typedef double Scalar;
typedef int    LocalOrdinal;
#ifdef HAVE_TEUCHOS_LONG_LONG_INT
typedef long long int GlobalOrdinal;
#else
typedef int GlobalOrdinal;
#endif
//
typedef Kokkos::DefaultNode::DefaultNodeType Node;
typedef Kokkos::DefaultKernels<Scalar,LocalOrdinal,Node>::SparseOps LocalMatOps;
//
#include "MueLu_UseShortNames.hpp"


template <typename Scalar, typename LocalOrdinal, typename GlobalOrdinal>
Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal> > TriDiag(const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal> > & map,
        const GlobalOrdinal nx, // note: nx unused
        const Scalar a, const Scalar b, const Scalar c)
{

  Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal> > mtx = Galeri::Xpetra::OperatorTraits<Xpetra::Map<LocalOrdinal, GlobalOrdinal>,Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal> >::Build(map, 3);

  LocalOrdinal NumMyElements = map->getNodeNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getNodeElementList();

  Teuchos::RCP<const Teuchos::Comm<int> > comm = map->getComm();

  GlobalOrdinal NumGlobalElements = map->getGlobalNumElements();

  GlobalOrdinal NumEntries;
  LocalOrdinal nnz=2;
  std::vector<Scalar> Values(nnz);
  std::vector<GlobalOrdinal> Indices(nnz);

  comm->barrier();

  Teuchos::RCP<Teuchos::Time> timer = rcp(new Teuchos::Time("TriDiag global insert"));
  timer->start(true);

  for (LocalOrdinal i = 0; i < NumMyElements; ++i)
    {
      if (MyGlobalElements[i] == 0)
        {
          // off-diagonal for first row
          Indices[0] = 1;
          NumEntries = 1;
          Values[0] = 0.0;//c; // dirichlet bc left (c)
        }
      else if (MyGlobalElements[i] == NumGlobalElements - 1)
        {
          // off-diagonal for last row
          Indices[0] = NumGlobalElements - 2;
          NumEntries = 1;
          Values[0] = 0.0; //Teuchos::ScalarTraits<Scalar>(0.0); // dirichlet bc right (b)
        }
      else
        {
          // off-diagonal for internal row
          Indices[0] = MyGlobalElements[i] - 1;
          Values[1] = b;
          Indices[1] = MyGlobalElements[i] + 1;
          Values[0] = c;
          NumEntries = 2;
        }

      // put the off-diagonal entries
      // Xpetra wants ArrayViews (sigh)
      Teuchos::ArrayView<Scalar> av(&Values[0],NumEntries);
      Teuchos::ArrayView<GlobalOrdinal> iv(&Indices[0],NumEntries);
      mtx->insertGlobalValues(MyGlobalElements[i], iv, av);

      // Put in the diagonal entry
      mtx->insertGlobalValues(MyGlobalElements[i],
                              Teuchos::tuple<GlobalOrdinal>(MyGlobalElements[i]),
                              Teuchos::tuple<Scalar>(a) );

    } //for (LocalOrdinal i = 0; i < NumMyElements; ++i)

    timer->stop();


  timer = rcp(new Teuchos::Time("TriDiag fillComplete"));
  timer->start(true);

  //mtx->fillComplete();

  timer->stop();

  return mtx;
} //TriDiag

int main(int argc, char *argv[]) {
  using Teuchos::RCP;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);

  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);
  *out << MueLu::MemUtils::PrintMemoryUsage() << std::endl;

  // Timing
  Teuchos::Time myTime("global");
  Teuchos::TimeMonitor M(myTime);
  Teuchos::Array<Teuchos::RCP<Teuchos::Time> > mtime;

  //out->setOutputToRootOnly(-1);
  //out->precision(12);

#ifndef HAVE_TEUCHOS_LONG_LONG_INT
  *out << "Warning: scaling test was not compiled with long long int support" << std::endl;
#endif

  /**********************************************************************************/
  /* SET TEST PARAMETERS                                                            */
  /**********************************************************************************/
  // Note: use --help to list available options.
  Teuchos::CommandLineProcessor clp(false);

  // Default is Laplace1D with nx = 8748.
  // It's a nice size for 1D and perfect aggregation. (6561=3^8)
  //Nice size for 1D and perfect aggregation on small numbers of processors. (8748=4*3^7)
  //Galeri::Xpetra::Parameters<GO> matrixParameters(clp, 8748); // manage parameters of the test case
  Xpetra::Parameters xpetraParameters(clp);             // manage parameters of xpetra

  // custom parameters
  // matrix parameters
  GO numGlobalElements = 1000;
  LO maxLevels = 5;
  LO its=10;
  std::string smooType="gs";
  Scalar smooDamping = 0.7;
  std::string transferOpType = "PG-AMG";
  int pauseForDebugger=0;
  int amgAsSolver=1;
  int amgAsPrecond=1;
  int sweeps=1;
  int maxCoarseSize=100;  //FIXME clp doesn't like long long int
  Scalar SADampingFactor=4./3;
  double tol = 1e-7;
  std::string aggOrdering = "natural";
  int minPerAgg=2;
  int maxNbrAlreadySelected=0;

  clp.setOption("numEle",&numGlobalElements,"number of elements (= dimension of problem - 1)");

  clp.setOption("maxLevels",&maxLevels,"maximum number of levels allowed");
  clp.setOption("its",&its,"number of multigrid cycles");
  clp.setOption("debug",&pauseForDebugger,"pause to attach debugger");
  clp.setOption("fixPoint",&amgAsSolver,"apply multigrid as solver");
  clp.setOption("precond",&amgAsPrecond,"apply multigrid as preconditioner");

  clp.setOption("transferOps",&transferOpType,"type of transfer operators (PA-AMG, SA-AMG, PG-AMG)");
  clp.setOption("saDamping",&SADampingFactor,"prolongator damping factor");

  clp.setOption("maxCoarseSize",&maxCoarseSize,"maximum #dofs in coarse operator");
  clp.setOption("tol",&tol,"stopping tolerance for Krylov method");
  clp.setOption("aggOrdering",&aggOrdering,"aggregation ordering strategy (natural,random,graph)");
  clp.setOption("minPerAgg",&minPerAgg,"minimum #DOFs per aggregate");
  clp.setOption("maxNbrSel",&maxNbrAlreadySelected,"maximum # of nbrs allowed to be in other aggregates");

  // level smoother settings
  clp.setOption("smooType",&smooType,"smoother type ('jacobi', 'gs', 'sgs', 'cheby')");
  clp.setOption("damping",&smooDamping,"damping factor for level smoothers");
  clp.setOption("sweeps",&sweeps,"sweeps to be used in level smoother (or Chebyshev degree)");

  switch (clp.parse(argc,argv)) {
  case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
  case Teuchos::CommandLineProcessor::PARSE_ERROR:
  case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
  case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  if (pauseForDebugger) {
    Utils::PauseForDebugger();
  }

  xpetraParameters.check();
  // TODO: check custom parameters
  std::transform(smooType.begin(), smooType.end(), smooType.begin(), ::tolower);
  Xpetra::UnderlyingLib lib = xpetraParameters.GetLib();

  if (comm->getRank() == 0) {
    std::cout << xpetraParameters;
    // TODO: print custom parameters
  }

  /**********************************************************************************/
  /* CREATE INITIAL MATRIX                                                          */
  /**********************************************************************************/
  mtime.push_back(M.getNewTimer("Matrix Build"));
  (mtime.back())->start();
  const RCP<const Map> map = MapFactory::Build(lib, numGlobalElements, 0, comm);
  RCP<Matrix> Op = TriDiag<SC,LO,GO>(map,0, 1.0, 0.0, -1.0);
  Op->fillComplete();
  mtime.back()->stop();

  // build RHS
  RCP<Vector> rhs = VectorFactory::Build(map);
  Scalar rhsval = (Scalar)1/numGlobalElements;
  rhs->putScalar(rhsval);

  Teuchos::ArrayRCP< Scalar > rhs_data = rhs->getDataNonConst(0);
  GlobalOrdinal NumGlobalElements = map->getGlobalNumElements();
  LocalOrdinal NumMyElements = map->getNodeNumElements();
  Teuchos::ArrayView<const GlobalOrdinal> MyGlobalElements = map->getNodeElementList();
  if (MyGlobalElements[0] == 0)  // left dirichlet bc
    rhs_data[0] = 0.0;
  if (MyGlobalElements[NumMyElements-1] == NumGlobalElements - 1) // right dirichlet bd
    rhs_data[NumMyElements-1] = 0.0;


  /**********************************************************************************/
  /* BUILD MG HIERARCHY                                                             */
  /**********************************************************************************/

  // build nullspace
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
  nullSpace->putScalar( (SC) 1.0);
  Teuchos::Array<ST::magnitudeType> norms(1);

  // MueLu setup
  mtime.push_back(M.getNewTimer("MueLu Setup"));
  mtime.back()->start(); // start time measurement
  RCP<MueLu::Hierarchy<SC,LO,GO,NO,LMO> > H = rcp ( new Hierarchy() );
  H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  H->SetMaxCoarseSize((GO) maxCoarseSize);;

  // build finest Level
  RCP<MueLu::Level> Finest = H->GetLevel();
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  Finest->Set("A",Op);
  Finest->Set("Nullspace",nullSpace);

  RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
  *out << "========================= Aggregate option summary  =========================" << std::endl;
  *out << "min DOFs per aggregate :                " << minPerAgg << std::endl;
  *out << "min # of root nbrs already aggregated : " << maxNbrAlreadySelected << std::endl;
  UCAggFact->SetMinNodesPerAggregate(minPerAgg);  //TODO should increase if run anything other than 1D
  UCAggFact->SetMaxNeighAlreadySelected(maxNbrAlreadySelected);
  std::transform(aggOrdering.begin(), aggOrdering.end(), aggOrdering.begin(), ::tolower);
  if (aggOrdering == "natural") {
       *out << "aggregate ordering :                    NATURAL" << std::endl;
       UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
  } else if (aggOrdering == "random") {
       *out << "aggregate ordering :                    RANDOM" << std::endl;
       UCAggFact->SetOrdering(MueLu::AggOptions::RANDOM);
  } else if (aggOrdering == "graph") {
       *out << "aggregate ordering :                    GRAPH" << std::endl;
       UCAggFact->SetOrdering(MueLu::AggOptions::GRAPH);
  } else {
    std::string msg = "main: bad aggregation option """ + aggOrdering + """.";
    throw(MueLu::Exceptions::RuntimeError(msg));
  }
  UCAggFact->SetPhase3AggCreation(0.5);
  *out << "=============================================================================" << std::endl;

  RCP<PFactory> Pfact = Teuchos::null;
  RCP<Factory> Rfact = Teuchos::null;

  if (transferOpType == "PA-AMG") {
    Pfact = rcp(new TentativePFactory(UCAggFact));
    Rfact = rcp( new TransPFactory(Pfact));
  }
  else if(transferOpType == "SA-AMG") {
    // build transfer operators
    RCP<TentativePFactory> TentPFact = rcp(new TentativePFactory(UCAggFact));
    Pfact = rcp( new SaPFactory(TentPFact) );
    Rfact = rcp( new TransPFactory(Pfact));
  }
  else if(transferOpType == "PG-AMG") {
    // build transfer operators
    RCP<TentativePFactory> TentPFact = rcp(new TentativePFactory(UCAggFact));
    Pfact = rcp( new PgPFactory(TentPFact) );
    Rfact = rcp( new GenericRFactory(Pfact));
  }

  RCP<RAPFactory>       Acfact = rcp( new RAPFactory(Pfact, Rfact) );
  Acfact->setVerbLevel(Teuchos::VERB_HIGH);

  *out << " after ACFactory " << std::endl;

  // build level smoothers

  RCP<SmootherPrototype> smooProto;
  std::string ifpackType;
  Teuchos::ParameterList ifpackList;
  ifpackList.set("relaxation: sweeps", (LO) sweeps);
  ifpackList.set("relaxation: damping factor", (SC) smooDamping); // 0.7
  if (smooType == "sgs") {
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Symmetric Gauss-Seidel");
  } else if (smooType == "gs") {
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Gauss-Seidel");
  } else if (smooType == "jacobi") {
    ifpackType = "RELAXATION";
    ifpackList.set("relaxation: type", "Jacobi");
  }
  else if (smooType == "cheby") {
    ifpackType = "CHEBYSHEV";
    ifpackList.set("chebyshev: degree", (LO) sweeps);
    ifpackList.set("chebyshev: ratio eigenvalue", (SC) 20);
    ifpackList.set("chebyshev: max eigenvalue", (double) -1.0);
    ifpackList.set("chebyshev: min eigenvalue", (double) 1.0);
    ifpackList.set("chebyshev: zero starting solution", true);
  }

  smooProto = rcp( new TrilinosSmoother(lib, ifpackType, ifpackList) );
  RCP<SmootherFactory> SmooFact;
  if (maxLevels > 1)
    SmooFact = rcp( new SmootherFactory(smooProto) );

  *out << " after SmootherFactory " << std::endl;

  Teuchos::ParameterList status;
  status = H->FullPopulate(*Pfact,*Rfact,*Acfact,*SmooFact,0,maxLevels);
  mtime.back()->stop(); // stop time measurement for MueLu setup
  H->SetCoarsestSolver(*SmooFact,MueLu::PRE);

  mtime.back()->stop();
  *out  << "======================\n Multigrid statistics \n======================" << std::endl;
  status.print(*out,Teuchos::ParameterList::PrintOptions().indent(2));

  /**********************************************************************************/
  /* SOLVE PROBLEM                                                                  */
  /**********************************************************************************/

  RCP<MultiVector> X = MultiVectorFactory::Build(map,1);

  // Use AMG directly as an iterative method
  {
    X->putScalar( (SC) 0.0);

    H->Iterate(*rhs,its,*X);

    //x->describe(*out,Teuchos::VERB_EXTREME);
  }

  // use MueLu as preconditioner within Belos
  if(amgAsPrecond && lib == Xpetra::UseTpetra)   // TODO: MueLu as preconditioner in Belos only works with Tpetra!
  {
#if defined(HAVE_MUELU_BELOS) && defined(HAVE_MUELU_TPETRA)
    X->putScalar( (SC) 0.0);

    int numrhs = 1;
    RCP<MultiVector> resid = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(map,numrhs);

    typedef ST::magnitudeType         MT;
    typedef Tpetra::MultiVector<SC>   MV;
    typedef Belos::OperatorT<MV>      OP;

    // Vectors
    RCP<MV> belosX     = MueLu::Utils<SC,LO,GO,NO,LMO>::MV2NonConstTpetraMV(X);
    RCP<MV> belosRHS   = MueLu::Utils<SC,LO,GO,NO,LMO>::MV2NonConstTpetraMV(rhs);
    RCP<MV> belosResid = MueLu::Utils<SC,LO,GO,NO,LMO>::MV2NonConstTpetraMV(resid);

    // construct Belos LinearProblem
    RCP<OP> belosOp      = Teuchos::rcp (new Belos::XpetraOp<SC,LO,GO,NO,LMO>(Op) );  // Xpetra::Op -> Belos::Op
    RCP<OP> belosPrec    = Teuchos::rcp (new Belos::MueLuOp<SC,LO,GO,NO,LMO>(H)); // Hierarchy  -> prec

    RCP<Belos::LinearProblem<double,MV,OP> > problem = Teuchos::rcp( new Belos::LinearProblem<double,MV,OP>(belosOp, belosX, belosRHS) );
    problem->setLeftPrec( belosPrec );

    bool set = problem->setProblem();
    if (set == false) {
      std::cout << std::endl << "ERROR: Belos::LinearProblem failed to set up correctly!" << std::endl;
      return EXIT_FAILURE;
    }

    // create an iterative solver manager
    int maxiters = 100;
    double tol = 1e-4;
    Teuchos::ParameterList belosList;
    belosList.set("Maximum Iterations", maxiters);
    belosList.set("Convergence Tolerance", tol);
    belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails);

    RCP< Belos::SolverManager<double,MV,OP> > solver = Teuchos::rcp( new Belos::BlockGmresSolMgr<double, MV, OP>(problem, rcp(&belosList,false)));

    // perform solve
    Belos::ReturnType ret;
    bool badRes = false;

    try{
      // Perform solve
      mtime.push_back(M.getNewTimer("Belos Solve"));
      mtime.back()->start();
      ret = solver->solve();
      mtime.back()->stop();

      // Get the number of iterations for this solve.
      int numIters = solver->getNumIters();
      *out << "Number of iterations performed for this solve: " << numIters << std::endl;

      // Compute actual residuals.
      std::vector<double> actual_resids( numrhs ); //TODO: double?
      std::vector<double> rhs_norm( numrhs );

      typedef Belos::OperatorTraits<SC,MV,OP>  OPT;
      typedef Belos::MultiVecTraits<SC,MV>     MVT;

      OPT::Apply( *belosOp, *belosX, *belosResid );
      MVT::MvAddMv( -1.0, *belosResid, 1.0, *belosRHS, *belosResid );
      MVT::MvNorm( *belosResid, actual_resids );
      MVT::MvNorm( *belosRHS, rhs_norm );
      *out<< "---------- Actual Residuals (normalized) ----------"<<std::endl<<std::endl;
      for ( int i=0; i<numrhs; i++) {
        double actRes = actual_resids[i]/rhs_norm[i];
        *out<<"Problem "<<i<<" : \t"<< actRes <<std::endl;
        if (actRes > tol) { badRes = true; }
      }
    } //try
    catch(...) {
      *out << std::endl << "ERROR:  Belos threw an error! " << std::endl;
    }

    // Check convergence
    if (ret!=Belos::Converged || badRes) {
      *out << std::endl << "ERROR:  Belos did not converge! " << std::endl;
    } else
      *out << std::endl << "SUCCESS:  Belos converged!" << std::endl;

#endif
  }


  // Final summaries - this eats memory like a hot dog eating contest
  // M.summarize();

  int ntimers=mtime.size();
  Teuchos::ArrayRCP<double> lTime(ntimers);
  Teuchos::ArrayRCP<double> gTime(ntimers);

  for(int i=0;i<ntimers;i++) lTime[i]=mtime[i]->totalElapsedTime();

  // Allreduce is my friend.
#ifdef HAVE_MPI
  MPI_Allreduce(&*lTime,&*gTime,ntimers,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
#else
  for(int i=0;i<ntimers;i++) gTime[i] = lTime[i];
#endif

  for(int i=0;i<ntimers;i++) *out<<mtime[i]->name()<<": \t"<<gTime[i]<<std::endl;

  *out << MueLu::MemUtils::PrintMemoryUsage() << std::endl;

  return EXIT_SUCCESS;
}



