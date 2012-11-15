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
 * visualizeAggregates.cpp
 *
 *  Created on: 19.10.2011
 *      Author: tobias
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

// Xpetra
#include <Xpetra_Map.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_VectorFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Parameters.hpp>

// MueLu
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Memory.hpp"
#include "MueLu_Hierarchy.hpp"
#include "MueLu_TentativePFactory.hpp"
#include "MueLu_SmootherFactory.hpp"
#include "MueLu_UCAggregationFactory.hpp"
#include "MueLu_PgPFactory.hpp"
#include "MueLu_GenericRFactory.hpp"
#include "MueLu_SaPFactory.hpp"
#include "MueLu_TransPFactory.hpp"
#include "MueLu_RAPFactory.hpp"
#include "MueLu_TrilinosSmoother.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_Aggregates.hpp"

#include "MueLu_UseDefaultTypes.hpp"
#include "MueLu_UseShortNames.hpp"

void ExportAggregates(const Teuchos::RCP<Level>& level, const MueLu::FactoryBase* aggFact,const Teuchos::RCP<const Teuchos::Comm<int> >& comm );

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
  Teuchos::TimeMonitor Monitor(myTime);


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
  Epetra_Map emap(10201,0,*Xpetra::toEpetra(comm));
  Epetra_CrsMatrix * ptrA = 0;
  Epetra_Vector * ptrf = 0;

  std::cout << "Reading matrix market file" << std::endl;
  EpetraExt::MatrixMarketFileToCrsMatrix("Condif2Mat.mat",emap,emap,emap,ptrA);
  EpetraExt::MatrixMarketFileToVector("Condif2Rhs.mat",emap,ptrf);
  RCP<Epetra_CrsMatrix> epA = Teuchos::rcp(ptrA);
  RCP<Epetra_Vector> epv = Teuchos::rcp(ptrf);

  // Epetra_CrsMatrix -> Xpetra::Matrix
  RCP<Xpetra::CrsMatrix<double, int, int> > exA = Teuchos::rcp(new Xpetra::EpetraCrsMatrix(epA));
  RCP<Xpetra::CrsMatrixWrap<double, int, int> > crsOp = Teuchos::rcp(new Xpetra::CrsMatrixWrap<double, int, int>(exA));
  RCP<Xpetra::Matrix<double, int, int> > Op = Teuchos::rcp_dynamic_cast<Xpetra::Matrix<double, int, int> >(crsOp);

  // Epetra_Vector -> Xpetra::Vector
  RCP<Xpetra::Vector<double,int,int> > xRhs = Teuchos::rcp(new Xpetra::EpetraVector(epv));

  // Epetra_Map -> Xpetra::Map
  const RCP< const Xpetra::Map<int, int> > map = Xpetra::toXpetra(emap);

  // build nullspace
  RCP<MultiVector> nullSpace = MultiVectorFactory::Build(map,1);
  nullSpace->putScalar( (SC) 1.0);

  RCP<MueLu::Hierarchy<SC,LO,GO,NO,LMO> > H = rcp ( new Hierarchy() );
  H->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  H->SetMaxCoarseSize((GO) maxCoarseSize);;

  // build finest Level
  RCP<MueLu::Level> Finest = H->GetLevel();
  Finest->setDefaultVerbLevel(Teuchos::VERB_HIGH);
  Finest->Set("A",Op);
  Finest->Set("Nullspace",nullSpace);

  RCP<UCAggregationFactory> UCAggFact = rcp(new UCAggregationFactory());
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
  Finest->Keep("Aggregates",UCAggFact.get());
  *out << "=============================================================================" << std::endl;

  // build transfer operators

  // build level smoothers
  RCP<SmootherPrototype> smooProto;
  std::string ifpackType;
  Teuchos::ParameterList ifpackList;
  ifpackList.set("relaxation: sweeps", (LO) 1);
  ifpackList.set("relaxation: damping factor", (SC) 0.9); // 0.7
  ifpackType = "RELAXATION";
  ifpackList.set("relaxation: type", "Gauss-Seidel");

  smooProto = Teuchos::rcp( new TrilinosSmoother(ifpackType, ifpackList) );
  RCP<SmootherFactory> SmooFact;
  if (maxLevels > 1)
    SmooFact = rcp( new SmootherFactory(smooProto) );

  RCP<FactoryManager> M = rcp(new FactoryManager());
  M->SetFactory("Aggregates", UCAggFact);
  M->SetFactory("Smoother", SmooFact);

  H->Setup(*M,0,maxLevels);

  // print out aggregation information
  for(LocalOrdinal l=0; l<H->GetNumLevels()-1;l++) {
    RCP<Level> level = H->GetLevel((int)l);
    ExportAggregates(level, UCAggFact.get(),comm);
  }

  return EXIT_SUCCESS;
}

void ExportAggregates(const Teuchos::RCP<Level>& level, const MueLu::FactoryBase* aggFact,const Teuchos::RCP<const Teuchos::Comm<int> >& comm ) {
  if(level->IsAvailable("Aggregates",aggFact) == false) {
    std::cout << "no Aggregates available." << std::endl;
    return;
  }
  if(level->IsAvailable("A") == false) {
    std::cout << "no matrix A available" << std::endl;
    return;
  }
  Teuchos::RCP<Aggregates> aggregates = level->Get< Teuchos::RCP<Aggregates> >("Aggregates",aggFact);
  Teuchos::RCP<Matrix> Op = level->Get<Teuchos::RCP<Matrix> >("A");

  Teuchos::RCP<LOVector>vertex2AggId_vector = aggregates->GetVertex2AggId();
  Teuchos::RCP<LOVector>procWinner_vector = aggregates->GetProcWinner();
  Teuchos::ArrayRCP<LO> vertex2AggId = aggregates->GetVertex2AggId()->getDataNonConst(0);
  Teuchos::ArrayRCP<LO> procWinner = aggregates->GetProcWinner()->getDataNonConst(0);

  // 1.) prepare for calculating global aggregate ids
  std::vector<GlobalOrdinal> numAggsGlobal(comm->getSize(),0);
  std::vector<GlobalOrdinal> numAggsLocal(comm->getSize(),0);
  std::vector<GlobalOrdinal> minGlobalAggId(comm->getSize(),0);
  numAggsLocal[comm->getRank()] = aggregates->GetNumAggregates();
  Teuchos::reduceAll(*comm,Teuchos::REDUCE_SUM, comm->getSize(),&numAggsLocal[0], &numAggsGlobal[0]);
  for(int i=1; i<Teuchos::as<int>(numAggsGlobal.size()); ++i) {
    numAggsGlobal[i] += numAggsGlobal[i-1];
    minGlobalAggId[i] = numAggsGlobal[i-1];
  }

  // 2.) build overlapping vector (global row id (column map) -> global AggId)
  Teuchos::RCP<MultiVector> OverlappingRowId2GlobalAggId = MultiVectorFactory::Build(Op->getColMap(),1);
  Teuchos::ArrayRCP<Scalar> OverlappingRowId2GlobalAggId_data = OverlappingRowId2GlobalAggId->getDataNonConst(0);

  // loop over local aggids
  for(LocalOrdinal lrow=0; lrow<vertex2AggId.size(); ++lrow)
  {
    LocalOrdinal localAggId = vertex2AggId[lrow];
    LocalOrdinal localAggIdProc = procWinner[lrow];
    GlobalOrdinal globalAggId = minGlobalAggId[localAggIdProc] + localAggId;
    OverlappingRowId2GlobalAggId_data[lrow] = globalAggId;
  }

  // 3.) communicate vector (global row id -> global agg id)
  Teuchos::RCP<const Import> importer = ImportFactory::Build(Op->getColMap(), Op->getRowMap());

  Teuchos::RCP<MultiVector> RowMap2GlobalAggId = MultiVectorFactory::Build(Op->getRowMap(),1);
  RowMap2GlobalAggId->doImport(*OverlappingRowId2GlobalAggId,*importer,Xpetra::INSERT);
  Teuchos::ArrayRCP<Scalar> RowMap2GlobalAggId_data = RowMap2GlobalAggId->getDataNonConst(0);

  // 4.) write to file
  std::stringstream filename;
  filename << "aggs_level" << level->GetLevelID() << "_proc" << comm->getRank() << ".out";
  std::ofstream fout(filename.str().c_str());
  for(size_t i=0; i<RowMap2GlobalAggId->getLocalLength(); i++) {
    fout << i << " " << Op->getRowMap()->getGlobalElement(i) << " " << RowMap2GlobalAggId_data[i] << std::endl;
  }
  fout.close();
}
