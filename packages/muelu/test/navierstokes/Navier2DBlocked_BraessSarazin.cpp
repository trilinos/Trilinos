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
/*
 * Navier2D_epetra.cpp
 *
 *  Created on: Mar 26, 2011
 *      Author: wiesner
 */

#include <unistd.h>
#include <iostream>
#include <fstream>

// Teuchos
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_ScalarTraits.hpp>

// Epetra
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_VectorIn.h>
#include <EpetraExt_MultiVectorIn.h>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_Parameters.hpp>
#include <Xpetra_MapExtractorFactory.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_StridedMapFactory.hpp>
#include <Xpetra_StridedMap.hpp>

// MueLu
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Memory.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_CoupledAggregationFactory.hpp"
#include "MueLu_UncoupledAggregationFactory.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_BlockedRAPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Aggregates.hpp"
#include "MueLu_CoalesceDropFactory.hpp"
#include "MueLu_PreDropFunctionConstVal.hpp"
#include "MueLu_NullspaceFactory.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_DirectSolver.hpp"
#include "MueLu_EpetraOperator.hpp"
#include "MueLu_SubBlockAFactory.hpp"
#include "MueLu_BlockedPFactory.hpp"
#include "MueLu_BlockedGaussSeidelSmoother.hpp"
#include "MueLu_SchurComplementFactory.hpp"
#include "MueLu_BraessSarazinSmoother.hpp"

#include "MueLu_CoarseMapFactory.hpp"

#include "MueLu_AmalgamationFactory.hpp"
#include "MueLu_AggregationExportFactory.hpp"

#include <Epetra_LinearProblem.h>
#include <AztecOO.h>

#include "Navier2D_Helpers.h"

/*!
 *  2d Navier Stokes example (for Epetra)
 *
 *  using block matrices
 */


int main(int argc, char *argv[]) {
#if defined(HAVE_MUELU_EPETRA) && defined(HAVE_MUELU_EPETRAEXT)
  typedef double Scalar;
  typedef int LocalOrdinal;
  typedef int GlobalOrdinal;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Xpetra::EpetraNode Node;
#include "MueLu_UseShortNames.hpp"

  using Teuchos::RCP;
  using Teuchos::rcp;
  using namespace MueLuTests;
  using namespace Teuchos;

  oblackholestream blackhole;
  GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  //
  RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  RCP<FancyOStream> out = fancyOStream(rcpFromRef(std::cout));
  out->setOutputToRootOnly(0);
  *out << MueLu::MemUtils::PrintMemoryUsage() << std::endl;

  // Timing
  Time myTime("global");
  TimeMonitor MM(myTime);

  // custom parameters
  LocalOrdinal maxLevels = 3;

  GlobalOrdinal maxCoarseSize=1; //FIXME clp doesn't like long long int

  int globalNumDofs = 1500;  // used for the maps
  int nDofsPerNode = 3;      // used for generating the fine level null-space

  // build strided maps
  // striding information: 2 velocity dofs and 1 pressure dof = 3 dofs per node
  std::vector<size_t> stridingInfo;
  stridingInfo.push_back(2);
  stridingInfo.push_back(1);

  /////////////////////////////////////// build strided maps
  // build strided maps:
  // xstridedfullmap: full map (velocity and pressure dof gids), continous
  // xstridedvelmap: only velocity dof gid maps (i.e. 0,1,3,4,6,7...)
  // xstridedpremap: only pressure dof gid maps (i.e. 2,5,8,...)
  Xpetra::UnderlyingLib lib = Xpetra::UseEpetra;
  RCP<StridedMap> xstridedfullmap = StridedMapFactory::Build(lib,globalNumDofs,0,stridingInfo,comm,-1);
  RCP<StridedMap> xstridedvelmap  = StridedMapFactory::Build(xstridedfullmap,0);
  RCP<StridedMap> xstridedpremap  = StridedMapFactory::Build(xstridedfullmap,1);

  /////////////////////////////////////// transform Xpetra::Map objects to Epetra
  // this is needed for AztecOO
  const RCP<const Epetra_Map> fullmap = rcpFromRef(Xpetra::toEpetra(*xstridedfullmap));
  RCP<const Epetra_Map>       velmap  = rcpFromRef(Xpetra::toEpetra(*xstridedvelmap));
  RCP<const Epetra_Map>       premap  = rcpFromRef(Xpetra::toEpetra(*xstridedpremap));

  /////////////////////////////////////// import problem matrix and RHS from files (-> Epetra)

  // read in problem
  Epetra_CrsMatrix * ptrA = 0;
  Epetra_Vector * ptrf = 0;
  Epetra_MultiVector* ptrNS = 0;

  *out << "Reading matrix market file" << std::endl;

  EpetraExt::MatrixMarketFileToCrsMatrix("A_re1000_5932.txt",*fullmap,*fullmap,*fullmap,ptrA);
  EpetraExt::MatrixMarketFileToVector("b_re1000_5932.txt",*fullmap,ptrf);

  RCP<Epetra_CrsMatrix> epA = rcp(ptrA);
  RCP<Epetra_Vector> epv = rcp(ptrf);
  RCP<Epetra_MultiVector> epNS = rcp(ptrNS);


  /////////////////////////////////////// split system into 2x2 block system

  *out << "Split matrix into 2x2 block matrix" << std::endl;

  // split fullA into A11,..., A22
  RCP<Epetra_CrsMatrix> A11;
  RCP<Epetra_CrsMatrix> A12;
  RCP<Epetra_CrsMatrix> A21;
  RCP<Epetra_CrsMatrix> A22;

  if(SplitMatrix2x2(epA,*velmap,*premap,A11,A12,A21,A22)==false)
    *out << "Problem with splitting matrix"<< std::endl;

  /////////////////////////////////////// transform Epetra objects to Xpetra (needed for MueLu)

  // build Xpetra objects from Epetra_CrsMatrix objects
  RCP<Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > xA11 = rcp(new Xpetra::EpetraCrsMatrixT<GlobalOrdinal,Node>(A11));
  RCP<Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > xA12 = rcp(new Xpetra::EpetraCrsMatrixT<GlobalOrdinal,Node>(A12));
  RCP<Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > xA21 = rcp(new Xpetra::EpetraCrsMatrixT<GlobalOrdinal,Node>(A21));
  RCP<Xpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > xA22 = rcp(new Xpetra::EpetraCrsMatrixT<GlobalOrdinal,Node>(A22));

  /////////////////////////////////////// generate MapExtractor object

  std::vector<RCP<const Xpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > > xmaps;
  xmaps.push_back(xstridedvelmap);
  xmaps.push_back(xstridedpremap);

  RCP<const Xpetra::MapExtractor<Scalar,LocalOrdinal,GlobalOrdinal,Node> > map_extractor = Xpetra::MapExtractorFactory<Scalar,LocalOrdinal,GlobalOrdinal,Node>::Build(xstridedfullmap,xmaps);

  /////////////////////////////////////// build blocked transfer operator
  // using the map extractor
  RCP<Xpetra::BlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> > bOp = rcp(new Xpetra::BlockedCrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node>(map_extractor,map_extractor,10));
  bOp->setMatrix(0,0,Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>(xA11)));
  bOp->setMatrix(0,1,Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>(xA12)));
  bOp->setMatrix(1,0,Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>(xA21)));
  bOp->setMatrix(1,1,Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar,LocalOrdinal,GlobalOrdinal,Node>(xA22)));

  bOp->fillComplete();

  //////////////////////////////////////////////////// create Hierarchy
  RCP<Hierarchy> H = rcp ( new Hierarchy() );
  H->setDefaultVerbLevel(VERB_HIGH);
  //H->setDefaultVerbLevel(VERB_NONE);
  H->SetMaxCoarseSize(maxCoarseSize);

  //////////////////////////////////////////////////////// finest Level
  RCP<MueLu::Level> Finest = H->GetLevel();
  Finest->setDefaultVerbLevel(VERB_HIGH);
  Finest->Set("A",rcp_dynamic_cast<Matrix>(bOp));

  /////////////////////////////////////////////// define subblocks of A
  // make A11 block and A22 block available as variable "A" generated
  // by A11Fact and A22Fact
  RCP<SubBlockAFactory> A11Fact = rcp(new SubBlockAFactory());
  A11Fact->SetFactory("A",MueLu::NoFactory::getRCP());
  A11Fact->SetParameter("block row",ParameterEntry(0));
  A11Fact->SetParameter("block col",ParameterEntry(0));
  RCP<SubBlockAFactory> A22Fact = rcp(new SubBlockAFactory());
  A22Fact->SetFactory("A",MueLu::NoFactory::getRCP());
  A22Fact->SetParameter("block row",ParameterEntry(1));
  A22Fact->SetParameter("block col",ParameterEntry(1));

  ////////////////////////////////////////// prepare null space for A11
  RCP<MultiVector> nullspace11 = MultiVectorFactory::Build(xstridedvelmap, 2);  // this is a 2D standard null space

  for (int i=0; i<nDofsPerNode-1; ++i) {
    ArrayRCP<Scalar> nsValues = nullspace11->getDataNonConst(i);
    int numBlocks = nsValues.size() / (nDofsPerNode - 1);
    for (int j=0; j< numBlocks; ++j) {
      nsValues[j*(nDofsPerNode - 1) + i] = 1.0;
    }
  }

  Finest->Set("Nullspace1",nullspace11);

  ///////////////////////////////////////// define CoalesceDropFactory and Aggregation for A11
  // set up amalgamation for A11. Note: we're using a default null space factory (null)
  RCP<AmalgamationFactory> amalgFact11 = rcp(new AmalgamationFactory());
  amalgFact11->SetFactory("A", A11Fact);

  amalgFact11->setDefaultVerbLevel(VERB_EXTREME);
  RCP<CoalesceDropFactory> dropFact11 = rcp(new CoalesceDropFactory());
  dropFact11->SetFactory("A", A11Fact);
  dropFact11->SetFactory("UnAmalgamationInfo", amalgFact11);
  dropFact11->setDefaultVerbLevel(VERB_EXTREME);
  //RCP<CoupledAggregationFactory> CoupledAggFact11 = rcp(new CoupledAggregationFactory());
  RCP<UncoupledAggregationFactory> CoupledAggFact11 = rcp(new UncoupledAggregationFactory());
  CoupledAggFact11->SetFactory("Graph", dropFact11);
  CoupledAggFact11->SetMinNodesPerAggregate(9);
  CoupledAggFact11->SetMaxNeighAlreadySelected(2);
  CoupledAggFact11->SetOrdering("natural");
  //CoupledAggFact11->SetPhase3AggCreation(0.5);

  ///////////////////////////////////////// define transfer ops for A11
#if 0
  // use PG-AMG
  RCP<PgPFactory> P11Fact = rcp(new PgPFactory());

  RCP<GenericRFactory> R11Fact = rcp(new GenericRFactory());
  RCP<NullspaceFactory> nspFact11 = rcp(new NullspaceFactory("Nullspace1",P11tentFact));

  RCP<NullspaceFactory> nspFact11 = rcp(new NullspaceFactory("Nullspace1"));

  RCP<CoarseMapFactory> coarseMapFact11 = rcp(new CoarseMapFactory());
  coarseMapFact11->setStridingData(stridingInfo);
  coarseMapFact11->setStridedBlockId(0);

  //////////////////////////////// define factory manager for (1,1) block
  RCP<FactoryManager> M11 = rcp(new FactoryManager());
  M11->SetFactory("A", A11Fact);
  M11->SetFactory("P", P11Fact);
  M11->SetFactory("R", R11Fact);
  M11->SetFactory("Aggregates", CoupledAggFact11);
  M11->SetFactory("UnAmalgamationInfo", amalgFact11);
  M11->SetFactory("Nullspace", nspFact11);
  // M11->SetFactory("Ptent", P11tentFact);
  M11->SetFactory("CoarseMap", coarseMapFact11);
#else
  RCP<TentativePFactory> P11Fact = rcp(new TentativePFactory());

  RCP<TransPFactory> R11Fact = rcp(new TransPFactory());

  RCP<NullspaceFactory> nspFact11 = rcp(new NullspaceFactory("Nullspace1"));
  nspFact11->SetFactory("Nullspace1",P11Fact);

  RCP<CoarseMapFactory> coarseMapFact11 = rcp(new CoarseMapFactory());
  coarseMapFact11->setStridingData(stridingInfo);
  coarseMapFact11->setStridedBlockId(0);

  //////////////////////////////// define factory manager for (1,1) block
  RCP<FactoryManager> M11 = rcp(new FactoryManager());
  M11->SetFactory("A", A11Fact);
  M11->SetFactory("P", P11Fact);
  M11->SetFactory("R", R11Fact);
  M11->SetFactory("Aggregates", CoupledAggFact11);
  M11->SetFactory("UnAmalgamationInfo", amalgFact11);
  M11->SetFactory("Nullspace", nspFact11);
  // M11->SetFactory("Ptent", P11Fact);
  M11->SetFactory("CoarseMap", coarseMapFact11);
  M11->SetFactory("Graph", dropFact11);
#endif
  M11->SetIgnoreUserData(true);               // always use data from factories defined in factory manager

  ////////////////////////////////////////// prepare null space for A22
  RCP<MultiVector> nullspace22 = MultiVectorFactory::Build(xstridedpremap, 1);  // this is a 2D standard null space
  ArrayRCP<Scalar> nsValues22 = nullspace22->getDataNonConst(0);
  for (int j=0; j< nsValues22.size(); ++j) {
    nsValues22[j] = 1.0;
  }

  Finest->Set("Nullspace2",nullspace22);

  ///////////////////////////////////////// define transfer ops for A22
#if 0
  // use PGAMG
  RCP<AmalgamationFactory> amalgFact22 = rcp(new AmalgamationFactory(A22Fact));
  RCP<TentativePFactory> P22tentFact = rcp(new TentativePFactory(CoupledAggFact11, amalgFact22));

  RCP<SaPFactory> P22Fact = rcp(new SaPFactory(P22tentFact));

  //RCP<GenericRFactory> R22Fact = rcp(new GenericRFactory(P22Fact));
  RCP<TransPFactory> R22Fact = rcp(new TransPFactory(P22Fact));

  RCP<NullspaceFactory> nspFact22 = rcp(new NullspaceFactory("Nullspace2",P22tentFact));
  RCP<CoarseMapFactory> coarseMapFact22 = rcp(new CoarseMapFactory(CoupledAggFact11, nspFact22));
  coarseMapFact22->setStridingData(stridingInfo);
  coarseMapFact22->setStridedBlockId(1);

  //////////////////////////////// define factory manager for (2,2) block
  RCP<FactoryManager> M22 = rcp(new FactoryManager());
  M22->SetFactory("A", A22Fact);
  M22->SetFactory("P", P22Fact);
  M22->SetFactory("R", R22Fact);
  M22->SetFactory("Aggregates", AggFact22);
  M22->SetFactory("Nullspace", nspFact22);
  M22->SetFactory("Ptent", P22tentFact);
  M22->SetFactory("CoarseMap", coarseMapFact22);
  M22->SetIgnoreUserData(true);               // always use data from factories defined in factory manager

#else
  // use TentativePFactory
  RCP<AmalgamationFactory> amalgFact22 = rcp(new AmalgamationFactory());
  RCP<TentativePFactory> P22Fact = rcp(new TentativePFactory()); // check me (fed with A22) wrong column GIDS!!!

  RCP<TransPFactory> R22Fact = rcp(new TransPFactory());

  RCP<NullspaceFactory> nspFact22 = rcp(new NullspaceFactory("Nullspace2"));
  nspFact22->SetFactory("Nullspace2", P22Fact);
  RCP<CoarseMapFactory> coarseMapFact22 = rcp(new CoarseMapFactory());
  coarseMapFact22->setStridingData(stridingInfo);
  coarseMapFact22->setStridedBlockId(1);

  //////////////////////////////// define factory manager for (2,2) block
  RCP<FactoryManager> M22 = rcp(new FactoryManager());
  M22->SetFactory("A", A22Fact);
  M22->SetFactory("P", P22Fact);
  M22->SetFactory("R", R22Fact);
  M22->SetFactory("Aggregates", CoupledAggFact11);
  M22->SetFactory("Nullspace", nspFact22);
  M22->SetFactory("UnAmalgamationInfo", amalgFact22); // TODO oops what about that? it was M11 before?
  M22->SetFactory("Ptent", P22Fact);
  M22->SetFactory("CoarseMap", coarseMapFact22);
  M22->SetIgnoreUserData(true);               // always use data from factories defined in factory manager
#endif

  /////////////////////////////////////////// define blocked transfer ops
  RCP<BlockedPFactory> PFact = rcp(new BlockedPFactory()); // use row map index base from bOp
  PFact->AddFactoryManager(M11);
  PFact->AddFactoryManager(M22);

  RCP<GenericRFactory> RFact = rcp(new GenericRFactory());
  RFact->SetFactory("P", PFact);

  RCP<Factory> AcFact = rcp(new BlockedRAPFactory());
  AcFact->SetFactory("P", PFact);
  AcFact->SetFactory("R", RFact);

  *out << "Creating Braess-Sarazin Smoother" << std::endl;

  //////////////////////////////////////////////////////////////////////
  // Smoothers

  //Another factory manager for braes sarazin smoother
  //Schur Complement Factory, using the factory to generate AcFact
  SC omega = 1.3;
  RCP<SchurComplementFactory> SFact = rcp(new SchurComplementFactory());
  SFact->SetParameter("omega", ParameterEntry(omega));
  SFact->SetFactory("A", MueLu::NoFactory::getRCP());

  //Smoother Factory, using SFact as a factory for A
  std::string ifpackSCType;
  ParameterList ifpackSCList;
  ifpackSCList.set("relaxation: sweeps", (LocalOrdinal) 3);
  ifpackSCList.set("relaxation: damping factor", (Scalar) 1.0);
  ifpackSCType = "RELAXATION";
  ifpackSCList.set("relaxation: type", "Gauss-Seidel");
  RCP<SmootherPrototype> smoProtoSC     = rcp( new TrilinosSmoother(ifpackSCType, ifpackSCList, 0) );
  smoProtoSC->SetFactory("A", SFact);
  RCP<SmootherFactory> SmooSCFact = rcp( new SmootherFactory(smoProtoSC) );

  RCP<BraessSarazinSmoother> smootherPrototype     = rcp( new BraessSarazinSmoother() );
  smootherPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(3));
  smootherPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(omega));

  RCP<SmootherFactory>   smootherFact          = rcp( new SmootherFactory(smootherPrototype) );

  RCP<BraessSarazinSmoother> coarseSolverPrototype = rcp( new BraessSarazinSmoother() );
  coarseSolverPrototype->SetParameter("Sweeps", Teuchos::ParameterEntry(3));
  coarseSolverPrototype->SetParameter("Damping factor", Teuchos::ParameterEntry(omega));

  RCP<SmootherFactory>   coarseSolverFact      = rcp( new SmootherFactory(coarseSolverPrototype, null) );

  RCP<FactoryManager> MB = rcp(new FactoryManager());
  MB->SetFactory("A",     SFact);
  MB->SetFactory("Smoother",    SmooSCFact);
  MB->SetIgnoreUserData(true);               // always use data from factories defined in factory manager
  smootherPrototype->AddFactoryManager(MB,0);
  coarseSolverPrototype->AddFactoryManager(MB,0);

  // main factory manager
  FactoryManager M;
  M.SetFactory("A",            AcFact);
  M.SetFactory("P",            PFact);
  M.SetFactory("R",            RFact);
  M.SetFactory("Smoother",     smootherFact); // TODO fix me
  M.SetFactory("PreSmoother",     smootherFact); // TODO fix me
  M.SetFactory("PostSmoother",     smootherFact); // TODO fix me
  M.SetFactory("CoarseSolver", coarseSolverFact);

  //////////////////////////////////// setup multigrid

  H->Setup(M,0,maxLevels);

  *out << std::endl;
  *out << "print content of multigrid levels:" << std::endl;

  Finest->print(*out);

  RCP<Level> coarseLevel = H->GetLevel(1);
  coarseLevel->print(*out);

  RCP<Level> coarseLevel2 = H->GetLevel(2);
  coarseLevel2->print(*out);

  RCP<MultiVector> xLsg = MultiVectorFactory::Build(xstridedfullmap,1);

  // Use AMG directly as an iterative method
#if 0
  {
    xLsg->putScalar( (SC) 0.0);

    // Epetra_Vector -> Xpetra::Vector
    RCP<Vector> xRhs = rcp(new Xpetra::EpetraVector(epv));

    // calculate initial (absolute) residual
    Array<ScalarTraits<SC>::magnitudeType> norms(1);
    xRhs->norm2(norms);
    *out << "||x_0|| = " << norms[0] << std::endl;

    // apply ten multigrid iterations
    H->Iterate(*xRhs,*xLsg,100);


    // calculate and print residual
    RCP<MultiVector> xTmp = MultiVectorFactory::Build(xstridedfullmap,1);
    bOp->apply(*xLsg,*xTmp,NO_TRANS,(SC)1.0,(SC)0.0);
    xRhs->update((SC)-1.0,*xTmp,(SC)1.0);
    xRhs->norm2(norms);
    *out << "||x|| = " << norms[0] << std::endl;
  }
#endif

  //
  // Solve Ax = b using AMG as a preconditioner in AztecOO
  //
  {
    RCP<Epetra_Vector> X = rcp(new Epetra_Vector(epv->Map()));
    X->PutScalar(0.0);
    Epetra_LinearProblem epetraProblem(epA.get(), X.get(), epv.get());

    AztecOO aztecSolver(epetraProblem);
    aztecSolver.SetAztecOption(AZ_solver, AZ_gmres);

    MueLu::EpetraOperator aztecPrec(H);
    aztecSolver.SetPrecOperator(&aztecPrec);

    int maxIts = 50;
    double tol = 1e-8;

    aztecSolver.Iterate(maxIts, tol);
  }

  return EXIT_SUCCESS;
#else
  std::cout << "Epetra (and/or EpetraExt) not enabled. Skip test." << std::endl;
  return EXIT_SUCCESS;
#endif // #if defined(HAVE_MUELU_SERIAL) && defined(HAVE_MUELU_EPETRA)
}
