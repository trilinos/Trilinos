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
 * Structure2D_epetra.cpp
 *
 *  Created on: Oct 24, 2011
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

// Epetra
#include <EpetraExt_CrsMatrixIn.h>
#include <EpetraExt_VectorIn.h>
#include <EpetraExt_MultiVectorIn.h>

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_Parameters.hpp>

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
#include "MueLu_Aggregates.hpp"
#include "MueLu_CoalesceFactory_def.hpp" // check me

// 
#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

/*!
 *  2d structural mechanics example for Epetra
 *
 *  (Nearly) Symmetric problem (except of Dirichlet boundaries) solved with AMG solver using a
 *  3 level multigrid with smoothed aggregation transfer operators.
 *
 */


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


#ifndef HAVE_TEUCHOS_LONG_LONG_INT
  *out << "Warning: scaling test was not compiled with long long int support" << std::endl;
#endif

  // custom parameters
  LO maxLevels = 4;

  int maxCoarseSize=1; //FIXME clp doesn't like long long int
  std::string aggOrdering = "natural";
  int minPerAgg=3;
  int maxNbrAlreadySelected=0;

  // read in problem
  Epetra_Map emap(3402,0,*Xpetra::toEpetra(comm));
  Epetra_CrsMatrix * ptrA = 0;
  Epetra_Vector * ptrf = 0;
  Epetra_MultiVector* ptrNS = 0;

  std::cout << "Reading matrix market file" << std::endl;
  EpetraExt::MatrixMarketFileToCrsMatrix("stru2d_A.txt",emap,emap,emap,ptrA);
  EpetraExt::MatrixMarketFileToVector("/home/tobias/trilinos/Trilinos_dev/ubuntu_openmpi/preCopyrightTrilinos/muelu/example/Structure/stru2d_b.txt",emap,ptrf);
  EpetraExt::MatrixMarketFileToMultiVector( "/home/tobias/trilinos/Trilinos_dev/ubuntu_openmpi/preCopyrightTrilinos/muelu/example/Structure/stru2d_ns.txt", emap, ptrNS);
  RCP<Epetra_CrsMatrix> epA = Teuchos::rcp(ptrA);
  RCP<Epetra_Vector> epv = Teuchos::rcp(ptrf);
  RCP<Epetra_MultiVector> epNS = Teuchos::rcp(ptrNS);

  // Epetra_CrsMatrix -> Xpetra::Matrix
  RCP<Xpetra::CrsMatrix<double, int, int> > exA = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(epA));
  RCP<Xpetra::CrsMatrixWrap<double, int, int> > crsOp = Teuchos::rcp(new Xpetra::CrsMatrixWrap<double, int, int>(exA));
  RCP<Xpetra::Matrix<double, int, int> > Op = Teuchos::rcp_dynamic_cast<Xpetra::Matrix<double, int, int> >(crsOp);

  // Epetra_Vector -> Xpetra::Vector
  RCP<Xpetra::Vector<double,int,int> > xRhs = Teuchos::rcp(new Xpetra::EpetraVector(epv));

  RCP<Xpetra::MultiVector<double,int,int> > xNS = Teuchos::rcp(new Xpetra::EpetraMultiVector(epNS));



  // Epetra_Map -> Xpetra::Map
  const RCP< const Xpetra::Map<int, int> > map = Xpetra::toXpetra(emap);

  RCP<MueLu::Hierarchy<SC,LO,GO,NO,LMO> > H = rcp ( new Hierarchy() );
  H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  H->SetMaxCoarseSize((GO) maxCoarseSize);;

  // build finest Level
  RCP<MueLu::Level> Finest = H->GetLevel();
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  Finest->Set("A",Op);
  Finest->Set("Nullspace",xNS);

  //RCP<CoalesceDropFactory> dropFact = rcp(new CoalesceDropFactory());
  RCP<CoalesceDropFactory> dropFact = rcp(new CoalesceDropFactory());
  RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory(dropFact));
  *out << "========================= Aggregate option summary  =========================" << std::endl;
  *out << "min DOFs per aggregate :                " << minPerAgg << std::endl;
  *out << "min # of root nbrs already aggregated : " << maxNbrAlreadySelected << std::endl;
  UCAggFact->SetMinNodesPerAggregate(minPerAgg); //TODO should increase if run anything other than 1D
  UCAggFact->SetMaxNeighAlreadySelected(maxNbrAlreadySelected);
  std::transform(aggOrdering.begin(), aggOrdering.end(), aggOrdering.begin(), ::tolower);
  if (aggOrdering == "natural") {
    *out << "aggregate ordering :                    NATURAL" << std::endl;
    UCAggFact->SetOrdering(MueLu::AggOptions::NATURAL);
  } else if (aggOrdering == "random") {
    *out << "aggregate ordering :                    RANDOM" << std::endl;
    UCAggFact->SetOrdering(MueLu::AggOptions::RANDOM);
  } else if (aggOrdering == "graph") {
    *out << "aggregate ordering :                    GRAPH" << std::endl;
    UCAggFact->SetOrdering(MueLu::AggOptions::GRAPH);
  } else {
    std::string msg = "main: bad aggregation option """ + aggOrdering + """.";
    throw(MueLu::Exceptions::RuntimeError(msg));
  }
  UCAggFact->SetPhase3AggCreation(0.5);
  *out << "=============================================================================" << std::endl;

  // build transfer operators
  RCP<NullspaceFactory> nspFact = rcp(new NullspaceFactory()); // make sure that we can keep nullspace!!!
  RCP<TentativePFactory> TentPFact = rcp(new TentativePFactory(UCAggFact,nspFact));
  //RCP<PgPFactory> Pfact = rcp( new PgPFactory(TentPFact) );
  //RCP<RFactory> Rfact  = rcp( new GenericRFactory(Pfact));
  RCP<SaPFactory> Pfact  = rcp( new SaPFactory(TentPFact) );
  RCP<RFactory>   Rfact  = rcp( new TransPFactory(Pfact) );
  RCP<RAPFactory> Acfact = rcp( new RAPFactory(Pfact, Rfact) );
  Acfact->setVerbLevel(Teuchos::VERB_HIGH);

  Finest->Keep("Aggregates",UCAggFact.get());
  Finest->Keep("Nullspace",nspFact.get());

  // build level smoothers
  RCP<SmootherPrototype> smooProto;
  std::string ifpackType;
  Teuchos::ParameterList ifpackList;
  ifpackList.set("relaxation: sweeps", (LO) 1);
  ifpackList.set("relaxation: damping factor", (SC) 0.9); // 0.7
  ifpackType = "RELAXATION";
  ifpackList.set("relaxation: type", "Gauss-Seidel");

  smooProto = Teuchos::rcp( new TrilinosSmoother(Xpetra::UseEpetra,ifpackType, ifpackList) );
  RCP<SmootherFactory> SmooFact;
  if (maxLevels > 1)
    SmooFact = rcp( new SmootherFactory(smooProto) );

  Teuchos::ParameterList status;
  status = H->FullPopulate(*Pfact,*Rfact,*Acfact,*SmooFact,0,maxLevels);

  H->SetCoarsestSolver(*SmooFact,MueLu::PRE);


  *out << "======================\n Multigrid statistics \n======================" << std::endl;
  status.print(*out,Teuchos::ParameterList::PrintOptions().indent(2));

  Finest->print(*out);

  RCP<Level> coarseLevel = H->GetLevel(1);
  coarseLevel->print(*out);

  RCP<Level> coarseLevel2 = H->GetLevel(2);
  coarseLevel2->print(*out);
  RCP<MultiVector> nsp2 = coarseLevel2->Get<RCP<MultiVector> >("Nullspace",nspFact.get());
  nsp2->describe(*out,Teuchos::VERB_EXTREME);

  RCP<Level> coarseLevel3 = H->GetLevel(3);
  coarseLevel3->print(*out);


  RCP<MultiVector> xLsg = MultiVectorFactory::Build(map,1);

  // Use AMG directly as an iterative method
  {
    xLsg->putScalar( (SC) 0.0);

    H->Iterate(*xRhs,10,*xLsg);

    //xLsg->describe(*out,Teuchos::VERB_EXTREME);
  }

  return EXIT_SUCCESS;
}
