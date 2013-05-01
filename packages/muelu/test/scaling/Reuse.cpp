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
#include <iostream>

#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_DefaultPlatform.hpp>

#include <Teuchos_Time.hpp>

// Galeri
#include <Galeri_XpetraParameters.hpp>
#include <Galeri_XpetraProblemFactory.hpp>
#include <Galeri_XpetraUtils.hpp>
#include <Galeri_XpetraMaps.hpp>
//

#include <MueLu.hpp>
#include <MueLu_Level.hpp>
#include <MueLu_BaseClass.hpp>
#include <MueLu_Exceptions.hpp>
#include <MueLu_ParameterListInterpreter.hpp> // TODO: move into MueLu.hpp

#include <MueLu_Utilities.hpp>

#include <MueLu_UseDefaultTypes.hpp>
#include <MueLu_UseShortNames.hpp>
#include <MueLu_MutuallyExclusiveTime.hpp>

#ifdef HAVE_MUELU_BELOS
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosXpetraAdapter.hpp>     // => This header defines Belos::XpetraOp
#include <BelosMueLuAdapter.hpp>      // => This header defines Belos::MueLuOp
#endif

using std::cout;
using std::endl;

int main(int argc, char *argv[]) {
  using Teuchos::RCP; // reference count pointers
  using Teuchos::rcp;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;

  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> sparse_matrix_type;
  typedef Tpetra::Map<LocalOrdinal, GlobalOrdinal, Node>                            map_type;
  typedef Teuchos::OrdinalTraits<Tpetra::global_size_t> TOT;
  typedef Teuchos::OrdinalTraits<LocalOrdinal> TOTLO;

  //
  // MPI initialization using Teuchos
  //

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
  RCP< const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int mypid = comm->getRank();

  RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

  //
  // Parameters
  //

  Teuchos::CommandLineProcessor clp(false); 
  Xpetra::Parameters             xpetraParameters(clp);                          // manage parameters of Xpetra

  std::string xmlFileName = "reuse.xml"; clp.setOption("xml",         &xmlFileName,  "read parameters from a file. Otherwise, this example uses by default 'reuse.xml'");
  std::string matrixPrefix = "jac"; clp.setOption("matrix",  &matrixPrefix,  "prefix for matrix file names.  Default = 'jac'");
  std::string rhsPrefix    = "rhs"; clp.setOption("rhs" ,    &matrixPrefix,  "prefix for rhs file names.  Default = 'rhs'");
  bool printTimings = true;         clp.setOption("timings", "notimings",  &printTimings, "print timings to screen");
  int first_matrix = 0;             clp.setOption("firstMatrix", &first_matrix, "first matrix in the sequence to use");
  int last_matrix = 1;              clp.setOption("lastMatrix",  &last_matrix,  "last matrix in the sequence to use");
  string do_reuse_str = "none";     clp.setOption("doReuse", &do_reuse_str, "if you want to try reuse");

  string matrixName;


  switch (clp.parse(argc,argv)) {
    case Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED:        return EXIT_SUCCESS; break;
    case Teuchos::CommandLineProcessor::PARSE_ERROR:
    case Teuchos::CommandLineProcessor::PARSE_UNRECOGNIZED_OPTION: return EXIT_FAILURE; break;
    case Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL:                               break;
  }

  int do_reuse=0;
  if(!strcmp(do_reuse_str.c_str(),"none")) do_reuse=0;
  else if(!strcmp(do_reuse_str.c_str(),"simple")) do_reuse=1;
  else if(!strcmp(do_reuse_str.c_str(),"fast")) do_reuse=2;
  else return EXIT_FAILURE;


  RCP<Matrix> Aprecond, Amatvec;
  RCP<MultiVector> rhs;

  RCP<TimeMonitor> globalTimeMonitor = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: S - Global Time")));
  RCP<TimeMonitor> tm, tm2;// = rcp (new TimeMonitor(*TimeMonitor::getNewTimer("ScalingTest: 1 - Matrix Build")));
  RCP<Time> timer;

  ParameterListInterpreter mueLuFactory(xmlFileName,*comm);

  // Operator and Multivector type that will be used with Belos
  typedef MultiVector          MV;
  typedef Belos::OperatorT<MV> OP;
 
  // Stats tracking
  int ArraySize = last_matrix-first_matrix+1;
  Array<Array<int>    > iteration_counts(ArraySize);
  Array<Array<double> > iteration_times(ArraySize);
  Array<Array<double> > setup_times(ArraySize);
  
  for(int i=0; i< ArraySize; i++) {
    iteration_counts[i].resize(ArraySize);
    iteration_times[i].resize(ArraySize);
    setup_times[i].resize(ArraySize);
  }



  for(int i=first_matrix; i <= last_matrix; i++) {
    char matrixFileName[80];
    char rhsFileName[80];
    char timerName[80];

    sprintf(matrixFileName,"%s%d.mm",matrixPrefix.c_str(),i);

    // Load the matrix 
    if(!mypid) std::cout<<"Loading matrix... "<<matrixFileName<<endl;
    Aprecond  = Utils::Read(string(matrixFileName), xpetraParameters.GetLib(), comm);
    Aprecond->SetFixedBlockSize(2);

    // Build the nullspace
    RCP<MultiVector> nullspace = MultiVectorFactory::Build(Aprecond->getRowMap(),2);
    nullspace->putScalar((Scalar)(0.0));

    Teuchos::ArrayRCP< Scalar > data0, data1;
    data0=nullspace->getDataNonConst(0); data1=nullspace->getDataNonConst(1);
    for(size_t k=0; k<Aprecond->getRowMap()->getNodeNumElements(); k++) {
      if( Aprecond->getRowMap()->getGlobalElement((LocalOrdinal)k) % 2 == 0)  data0[k]=1.0;
      else data1[k]=1.0;
    }
      
    // Build the preconditioner
    if(!mypid) std::cout<<"Building preconditioner... "<<matrixFileName<<endl;
    sprintf(timerName,"Reuse: Preconditioner Setup i=%d",i);
    tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer(timerName)));

    RCP<Hierarchy> H = mueLuFactory.CreateHierarchy();
    H->SetDefaultVerbLevel(MueLu::Extreme);
    H->GetLevel(0)->Set("A", Aprecond);
    H->GetLevel(0)->Set("Nullspace", nullspace);
    H->IsPreconditioner(true);

    if(do_reuse==2) {
      // Flag some things as keepers.
      H->GetLevel(0)->Keep("AP Pattern",mueLuFactory.GetFactoryManager(0)->GetFactory("A").get());
      //      H->GetLevel(0)->Keep("RAP Pattern",mueLuFactory.GetFactoryManager(0)->GetFactory("A").get());
    }

    mueLuFactory.SetupHierarchy(*H);    


    Teuchos::RCP<OP> belosPrec = Teuchos::rcp(new Belos::MueLuOp<SC, LO, GO, NO, LMO>(H));  // Turns a MueLu::Hierarchy object into a Belos operator
    tm=Teuchos::null;

    // Loop over all future matrices
    for(int j=i; j <= last_matrix; j++) {
      sprintf(matrixFileName,"%s%d.mm",matrixPrefix.c_str(),j);
      sprintf(rhsFileName,"%s%d.mm",rhsPrefix.c_str(),j);

      if(j != i) {
	// Load the matrix 
	if(!mypid) std::cout<<"-Loading matrix "<<matrixFileName<<endl;
	Amatvec  = Utils::Read(string(matrixFileName), xpetraParameters.GetLib(), comm);
      }
      else 
	Amatvec = Aprecond;
      Amatvec->SetFixedBlockSize(2);

      // Preconditioner update
      sprintf(timerName,"Reuse: Preconditioner Update i=%d",i);
      tm = rcp (new TimeMonitor(*TimeMonitor::getNewTimer(timerName)));
      // No-op at present

      if(do_reuse==0 && j!=i) {
	// No reuse: Do a full recompute
	H->GetLevel(0)->Set("A", Amatvec);
	mueLuFactory.SetupHierarchy(*H);
      }
      else if(do_reuse==2 && j!=i) {
	// "Fast" reuse
	// NTS: At the moment, this is equivalent to a full recompute
	H->GetLevel(0)->Set("A", Amatvec);
	mueLuFactory.SetupHierarchy(*H);
      }

      tm = Teuchos::null;

      // Load the RHS
      if(!mypid) std::cout<<"-Loading rhs "<<rhsFileName<<endl;
      rhs      = Utils::Read(string(rhsFileName), Amatvec->getRowMap());

      // Create an LHS
      RCP<Vector> X = VectorFactory::Build(Amatvec->getRowMap()); 
      X->putScalar(0.0);
  
      // Define Operator and Preconditioner
      Teuchos::RCP<OP> belosOp   = Teuchos::rcp(new Belos::XpetraOp<SC, LO, GO, NO, LMO>(Amatvec)); // Turns a Xpetra::Matrix object into a Belos operator
  
      // Construct a Belos LinearProblem object
      RCP< Belos::LinearProblem<SC, MV, OP> > belosProblem = rcp(new Belos::LinearProblem<SC, MV, OP>(belosOp, X, rhs));
      belosProblem->setLeftPrec(belosPrec);
      belosProblem->setProblem();


      // Belos parameter list
      int maxIts = 100;
      double tol = 1e-12;
      Teuchos::ParameterList belosList;
      belosList.set("Maximum Iterations",    maxIts); // Maximum number of iterations allowed
      belosList.set("Convergence Tolerance", tol);    // Relative convergence tolerance requested
      //belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::TimingDetails + Belos::StatusTestDetails);
      belosList.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
      belosList.set("Output Frequency",1);
      belosList.set("Output Style",Belos::Brief);
      
      // Create an iterative solver manager
      RCP< Belos::SolverManager<SC, MV, OP> > solver = rcp(new Belos::BlockCGSolMgr<SC, MV, OP>(belosProblem, rcp(&belosList, false)));

      // Perform solve
      sprintf(timerName,"Reuse: Solve i=%d j=%d",i,j);
      timer = TimeMonitor::getNewTimer(timerName);
      timer->start();
      Belos::ReturnType ret=Belos::Unconverged;
      ret = solver->solve();
      double my_time=timer->stop();	
      timer=Teuchos::null;

      // Get the number of iterations for this solve.
      if (comm->getRank() == 0)
	std::cout << "Number of iterations performed for this solve: " << solver->getNumIters() << std::endl;
      iteration_counts[i-first_matrix][j-first_matrix] = solver->getNumIters();
      iteration_times[i-first_matrix][j-first_matrix]  = my_time;

      // Check convergence
      if (ret != Belos::Converged) {
	if (!mypid) std::cout << std::endl << "ERROR:  Belos did not converge! " << std::endl;
      } else {
	if (!mypid) std::cout << std::endl << "SUCCESS:  Belos converged!" << std::endl;
      }

    }//end j
  }// end i

 

  globalTimeMonitor = Teuchos::null;
  
  if (printTimings)
    TimeMonitor::summarize(comm.ptr(), std::cout, false, true, false, Teuchos::Union);
  
 
  if(!mypid) {
    printf("************************* Iteration Counts ***********************\n");
    for(int i=0; i< ArraySize; i++) {
      for(int j=0; j< ArraySize; j++) 
	printf("%3d ",iteration_counts[i][j]);
      printf("\n");
    }
    
    printf("************************* Iteration Times ***********************\n");
    for(int i=0; i< ArraySize; i++) {
      for(int j=0; j< ArraySize; j++) 
	printf("%10.2f ",iteration_times[i][j]);
      printf("\n");
    }

   printf("************************* Setup Times ***********************\n");
    for(int i=0; i< ArraySize; i++) {
      for(int j=0; j< ArraySize; j++) 
	printf("%10.2f ",setup_times[i][j]);
      printf("\n");
    }


  }

//MueLu::MutuallyExclusiveTime<MueLu::BaseClass>::PrintParentChildPairs();


} //main
