//@HEADER
// ************************************************************************
//
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
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
// Questions? Contact Alexander Heinlein (alexander.heinlein@uni-koeln.de)
//
// ************************************************************************
//@HEADER

//#define Tpetra_issue_1752

#include <mpi.h>
#include <Epetra_MpiComm.h>

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>

#include <Galeri_Maps.h>
#include <Galeri_CrsMatrices.h>

#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_EpetraCrsMatrix.hpp>
#include <Xpetra_MapFactory.hpp>

#include <BelosOperatorT.hpp>
#include <BelosXpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>
//#include <BelosPseudoBlockGmresSolMgr.hpp>

#include <FROSch_TwoLevelBlockPreconditioner_def.hpp>

//#include "Tools/FROSch_Tools_def.hpp"

typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef int GO;
typedef KokkosClassic::DefaultNode::DefaultNodeType EpetraNode;
typedef EpetraNode NO;

typedef Xpetra::BlockedCrsMatrix<SC, LO, GO, NO> BlockedCrsMatrixClass;
typedef Xpetra::BlockedMap<LO, GO, NO> BlockedMapClass;
typedef Xpetra::Map<LO, GO, NO> MapClass;
typedef Xpetra::MultiVector< SC, LO, GO, NO > MultiVectorClass;

using namespace std;
using namespace Teuchos;
using namespace Xpetra;
using namespace FROSch;
using namespace Belos;

int main(int argc, char *argv[])
{
	MPI_Init(&argc,&argv);

	{

		Epetra_MpiComm CommWorld(MPI_COMM_WORLD);

		CommandLineProcessor My_CLP;

		int M = 4;
		My_CLP.setOption("M",&M,"H / h.");
		int Dimension = 2;
		My_CLP.setOption("Dim",&Dimension,"Dimension.");
		int Overlap = 1;
		My_CLP.setOption("Overlap",&Overlap,"Overlap.");
		int DofsPerNode = 1;
		My_CLP.setOption("DPN",&DofsPerNode,"Dofs per node.");
		int DOFOrdering = 0;
		My_CLP.setOption("Ordering",&DOFOrdering,"Dofs ordering (NodeWise=0,DimensionWise=1,Custom=2).");
		string xmlFile = "GDSW.xml";
		My_CLP.setOption("List",&xmlFile,"File name of the parameter list.");

		My_CLP.recogniseAllOptions(true);
		My_CLP.throwExceptions(false);
		CommandLineProcessor::EParseCommandLineReturn parseReturn = My_CLP.parse(argc,argv);
		if(parseReturn == CommandLineProcessor::PARSE_HELP_PRINTED) {
			MPI_Finalize();
			return 0;
		}

		int N;
		MPI_Comm COMM;
		int color=1;
		//bool onFirstLevelComm=false;
		if (Dimension == 2) {
			N = (int) (pow(CommWorld.NumProc(),1/2.) + 100*numeric_limits<double>::epsilon()); // 1/H
			if (CommWorld.MyPID()<N*N) {
				color=0;
			}
		} else if (Dimension == 3) {
			N = (int) (pow(CommWorld.NumProc(),1/3.) + 100*numeric_limits<double>::epsilon()); // 1/H
			if (CommWorld.MyPID()<N*N*N) {
				color=0;
			}
		} else {
			assert(0!=0);
		}

		MPI_Comm_split(CommWorld.Comm(),color,CommWorld.MyPID(),&COMM);
		RCP<Epetra_MpiComm> Comm(new Epetra_MpiComm(COMM));
		RCP<const Teuchos::Comm<int> > TeuchosComm = rcp(new MpiComm<int> (COMM));
		if (color==0) {

			RCP<ParameterList> parameterList = getParametersFromXmlFile(xmlFile);
            
//            Teuchos::RCP<Xpetra::Map<LO,GO,NO> > map = Xpetra::MapFactory<LO,GO,NO>::Build(UseTpetra,-1,1,0,TeuchosComm);
//            
//            Teuchos::Array<GO> small(0);
//            Teuchos::Array<GO> larger(0);
//            
//            small.push_back(TeuchosComm->getRank());
//            small.push_back((TeuchosComm->getRank()+1)%4);
//            
//            for (unsigned i=0; i<TeuchosComm->getSize(); i++) {
//                larger.push_back(i);
//            }
//            
//            Teuchos::RCP<Xpetra::Map<LO,GO,NO> > mapOlSmall = Xpetra::MapFactory<LO,GO,NO>::Build(UseTpetra,-1,small(),0,TeuchosComm);
//            Teuchos::RCP<Xpetra::Map<LO,GO,NO> > mapOlLarge = Xpetra::MapFactory<LO,GO,NO>::Build(UseTpetra,-1,larger(),0,TeuchosComm);
//            
//            SchwarzOperator<SC,LO,GO,NO>::MultiVectorPtr xOverlapSmall = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(mapOlSmall,1);
//            SchwarzOperator<SC,LO,GO,NO>::MultiVectorPtr xOverlapLarge = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(mapOlLarge,1);
//            
//            SchwarzOperator<SC,LO,GO,NO>::MultiVectorPtr x = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(map,1);
//            
//            xOverlapSmall->putScalar(1.);
//            xOverlapLarge->putScalar(1.);
//            
//            SchwarzOperator<SC,LO,GO,NO>::ImporterPtr Scatter = Xpetra::ImportFactory<LO,GO,NO>::Build(map,mapOlLarge);
//
//            SchwarzOperator<SC,LO,GO,NO>::ImporterPtr Scatter2 = Xpetra::ImportFactory<LO,GO,NO>::Build(map,mapOlSmall);
//            RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(std::cout));
//            xOverlapLarge->describe(*fancy,Teuchos::VERB_EXTREME);
//            mapOlSmall->describe(*fancy,Teuchos::VERB_EXTREME);
//            
//            x->putScalar(0.);
//            x->doExport(*xOverlapLarge,*Scatter2,Xpetra::ADD);
//            x->describe(*fancy,Teuchos::VERB_EXTREME);
//
//            x->putScalar(0.);
//            x->doExport(*xOverlapLarge,*Scatter,Xpetra::ADD);
//            x->describe(*fancy,Teuchos::VERB_EXTREME);
//
//            SchwarzOperator<SC,LO,GO,NO>::ImporterPtr ScatterSmallLarge = Xpetra::ImportFactory<LO,GO,NO>::Build(mapOlSmall,mapOlLarge);
//            xOverlapSmall->putScalar(0.);
//
//            SchwarzOperator<SC,LO,GO,NO>::MultiVectorPtr xOverlapSmallTmp = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(mapOlSmall,1);
//            xOverlapSmallTmp->doExport(*xOverlapLarge,*ScatterSmallLarge,Xpetra::INSERT);
//            xOverlapSmallTmp->describe(*fancy,Teuchos::VERB_EXTREME);
//            
//            x->putScalar(0.);
//            x->doExport(*xOverlapLarge,*Scatter2,Xpetra::ADD);
//            x->describe(*fancy,Teuchos::VERB_EXTREME);

            
			if (Comm->MyPID()==0) {
				cout << "--------------------------------------------------------------------------------\nPARAMETERS:" << endl;
				parameterList->print(cout);
				cout << "--------------------------------------------------------------------------------\n\n";
			}

            
			if (Comm->MyPID()==0) cout << "ASSEMBLY...";

			ParameterList GalerList;
			GalerList.set("nx", N*M);
			GalerList.set("ny", N*M);
			GalerList.set("nz", N*M);
			GalerList.set("mx", N);
			GalerList.set("my", N);
			GalerList.set("mz", N);

			RCP<Epetra_Map> UniqueMapEpetra;
			RCP<Epetra_CrsMatrix> KEpetra1;
            RCP<Epetra_CrsMatrix> KEpetra2;
            
			if (Dimension==2) {
				UniqueMapEpetra.reset(Galeri::CreateMap("Cartesian2D", *Comm, GalerList));
				KEpetra1.reset(Galeri::CreateCrsMatrix("Laplace2D", UniqueMapEpetra.get(), GalerList));
                KEpetra2.reset(Galeri::CreateCrsMatrix("Laplace2D", UniqueMapEpetra.get(), GalerList));
			} else if (Dimension==3) {
//				UniqueMapEpetra.reset(Galeri::CreateMap("Cartesian3D", *Comm, GalerList));
//				KEpetra.reset(Galeri::CreateCrsMatrix("Laplace3D", UniqueMapEpetra.get(), GalerList));
			}

			RCP<Map<LO,GO,NO> > UniqueMap1;
			RCP<Map<LO,GO,NO> > UniqueMap2;
			RCP<Matrix<SC,LO,GO,NO> > K1;
			RCP<Matrix<SC,LO,GO,NO> > K2;
			DofOrdering Ord;

			if (DOFOrdering == 0) {
				Ord = NodeWise;
				Array<GO> uniqueMapArray(DofsPerNode*UniqueMapEpetra->NumMyElements());
				for (LO i=0; i<UniqueMapEpetra->NumMyElements(); i++) {
					for (LO j=0; j<DofsPerNode; j++) {
						uniqueMapArray[DofsPerNode*i+j] = DofsPerNode*UniqueMapEpetra->GID(i)+j;
					}
				}

				UniqueMap1 = MapFactory<LO,GO,NO>::Build(UseTpetra,-1,uniqueMapArray(),0,TeuchosComm);
				K1 = MatrixFactory<SC,LO,GO,NO>::Build(UniqueMap1,KEpetra1->MaxNumEntries());
				K2 = MatrixFactory<SC,LO,GO,NO>::Build(UniqueMap1,KEpetra2->MaxNumEntries());
				for (LO i=0; i<UniqueMapEpetra->NumMyElements(); i++) {
					LO numEntries;
					GO* indices;
					SC* values;
					KEpetra1->ExtractMyRowView(i,numEntries,values,indices);

					for (LO j=0; j<DofsPerNode; j++) {
						Array<GO> indicesArray(numEntries);
						ArrayView<SC> valuesArrayView(values,numEntries);
						for (LO k=0; k<numEntries; k++) {
							indicesArray[k] = DofsPerNode*KEpetra1->ColMap().GID(indices[k])+j;
						}
						K1->insertGlobalValues(DofsPerNode*KEpetra1->RowMap().GID(i)+j,indicesArray(),valuesArrayView);
                        K2->insertGlobalValues(DofsPerNode*KEpetra2->RowMap().GID(i)+j,indicesArray(),valuesArrayView);
					}
				}
				K1->fillComplete();
				K2->fillComplete();
			} else if (DOFOrdering == 1) {
                assert(0!=0);
//				Ord = DimensionWise;
//				Array<GO> uniqueMapArray(DofsPerNode*UniqueMapEpetra->NumMyElements());
//				for (LO i=0; i<UniqueMapEpetra->NumMyElements(); i++) {
//					for (LO j=0; j<DofsPerNode; j++) {
//						uniqueMapArray[i+UniqueMapEpetra->NumMyElements()*j] = UniqueMapEpetra->GID(i)+(UniqueMapEpetra->MaxAllGID()+1)*j;
//					}
//				}
//
//				UniqueMap = MapFactory<LO,GO,NO>::Build(UseTpetra,-1,uniqueMapArray(),0,TeuchosComm);
//				K = MatrixFactory<SC,LO,GO,NO>::Build(UniqueMap,KEpetra->MaxNumEntries());
//				for (LO i=0; i<UniqueMapEpetra->NumMyElements(); i++) {
//					LO numEntries;
//					GO* indices;
//					SC* values;
//					KEpetra->ExtractMyRowView(i,numEntries,values,indices);
//
//					for (LO j=0; j<DofsPerNode; j++) {
//						Array<GO> indicesArray(numEntries);
//						ArrayView<SC> valuesArrayView(values,numEntries);
//						for (LO k=0; k<numEntries; k++) {
//							indicesArray[k] = KEpetra->ColMap().GID(indices[k])+(KEpetra->ColMap().MaxAllGID()+1)*j;
//						}
//						K->insertGlobalValues(UniqueMapEpetra->GID(i)+(UniqueMapEpetra->MaxAllGID()+1)*j,indicesArray(),valuesArrayView);
//					}
//				}
//				K->fillComplete();
			} else if (DOFOrdering == 2) {
				assert(0!=0); // TODO: Andere Sortierung implementieren
			} else {
				assert(0!=0);
			}
           
            
            
            RCP<const Teuchos::Comm<int> > CommX = K1->getRowMap()->getComm();
            UnderlyingLib lib = UniqueMap1->lib();
            
            if (Comm->MyPID()==0) cout << "done" << endl << "CONSTRUCTING BLOCKED MATRIX...";
            
            std::vector<RCP<const MapClass> > mapVector(2,Teuchos::null);
            
            mapVector.at(0) = K1->getRowMap();
            GO maxGIDPreviousMap = K1->getRowMap()->getMaxAllGlobalIndex();
            GO numGlobalElements = maxGIDPreviousMap + 1;
            Array< GO > elementList(0);
            Array< GO > elementListFull( K1->getRowMap()->getNodeElementList() );
            for (UN i=0; i<K2->getRowMap()->getNodeNumElements(); i++) {
                elementList.push_back(K2->getRowMap()->getGlobalElement(i) + numGlobalElements);
                elementListFull.push_back(K2->getRowMap()->getGlobalElement(i) + numGlobalElements);
            }
            mapVector.at(1) = MapFactory<LO, GO, NO>::Build(lib, numGlobalElements, elementList, 0, CommX);

            RCP<const MapClass > mapFull = MapFactory<LO, GO, NO>::Build(lib, 2*numGlobalElements, elementListFull, 0, CommX);
            RCP< MapClass > mapFullNonConst = MapFactory<LO, GO, NO>::Build(lib, 2*numGlobalElements, elementListFull, 0, CommX);

            RCP<Matrix<SC,LO,GO,NO> > KFull = MatrixFactory<SC,LO,GO,NO>::Build(mapFullNonConst,KEpetra1->MaxNumEntries());
            for (LO i=0; i<UniqueMapEpetra->NumMyElements(); i++) {
                LO numEntries;
                GO* indices;
                SC* values;
                KEpetra1->ExtractMyRowView(i,numEntries,values,indices);
                
                for (LO j=0; j<DofsPerNode; j++) {
                    Array<GO> indicesArray(numEntries);
                    ArrayView<SC> valuesArrayView(values,numEntries);
                    for (LO k=0; k<numEntries; k++) {
                        indicesArray[k] = DofsPerNode*KEpetra1->ColMap().GID(indices[k])+j;
                    }
                    KFull->insertGlobalValues(DofsPerNode*KEpetra1->RowMap().GID(i)+j,indicesArray(),valuesArrayView);
                }
                for (LO j=0; j<DofsPerNode; j++) {
                    Array<GO> indicesArray(numEntries);
                    ArrayView<SC> valuesArrayView(values,numEntries);
                    for (LO k=0; k<numEntries; k++) {
                        indicesArray[k] = DofsPerNode*KEpetra2->ColMap().GID(indices[k])+j + UniqueMap1->getMaxAllGlobalIndex() + 1;
                    }
                    KFull->insertGlobalValues(UniqueMap1->getMaxAllGlobalIndex() + 1+DofsPerNode*KEpetra1->RowMap().GID(i)+j,indicesArray(),valuesArrayView);
                }
            }
            KFull->fillComplete();
            RCP<const BlockedMapClass > mapBlocked(new BlockedMapClass(mapFull,mapVector));
            
            RCP<BlockedCrsMatrixClass> blockedMat(new BlockedCrsMatrixClass( mapBlocked, mapBlocked, 10 ));

            RCP<Matrix<SC,LO,GO,NO> > K3;
            RCP<Matrix<SC,LO,GO,NO> > K4;
            
            blockedMat->setMatrix (0, 0, K1);
            blockedMat->setMatrix (1, 1, K2);
//            blockedMat->setMatrix (0, 1, K3);
//            blockedMat->setMatrix (1, 0, K4);
            blockedMat->fillComplete();


            RCP<  MultiVectorClass > xFull = MultiVectorFactory<SC, LO, GO, NO>::Build(mapFull , 1);
            RCP<  MultiVectorClass > yFull = MultiVectorFactory<SC, LO, GO, NO>::Build(mapFull , 1);
//            xFull->putScalar(1.);
            yFull->putScalar(1.);
//            blockedMat->apply(*xFull,*yFull);
//            yFull->describe(*fancy,Teuchos::VERB_EXTREME);
            
//            Teuchos::RCP< Matrix<SC,LO,GO,NO> > system = blockedMat->Merge();
//            system->describe(*fancy,Teuchos::VERB_HIGH);
//            system->apply(*xFull,*yFull);
//            yFull->describe(*fancy,Teuchos::VERB_EXTREME);
            
			if (Comm->MyPID()==0) cout << "done" << endl << "CONSTRUCTING PRECONDITIONER...";

//			RCP<OneLevelPreconditioner<SC,LO,GO,NO> > OneLevelPrec(new OneLevelPreconditioner<SC,LO,GO,NO>(KFull , sublist(parameterList,"TwoLevelPreconditioner")));
//            OneLevelPrec->initialize(Overlap,mapFullNonConst);
            
            if (false) {
                RCP<OneLevelPreconditioner<SC,LO,GO,NO> > OneLevelPrec(new OneLevelPreconditioner<SC,LO,GO,NO>(K1 , sublist(parameterList,"TwoLevelPreconditioner")));
                if (Comm->MyPID()==0) cout << "INITIALIZE...";
                OneLevelPrec->initialize(Overlap);
                if (Comm->MyPID()==0) cout << "COMPUTE...";
                OneLevelPrec->compute();

                if (Comm->MyPID()==0) cout << "done" << endl << "SOLVING EQUATION SYSTEM...";
                
                RCP<MultiVector<SC,LO,GO,NO> > xSolution = MultiVectorFactory<SC,LO,GO,NO>::Build(UniqueMap1,1);
                RCP<MultiVector<SC,LO,GO,NO> > xRightHandSide = MultiVectorFactory<SC,LO,GO,NO>::Build(UniqueMap1,1);
                
                xSolution->putScalar(0.0);
                xRightHandSide->putScalar(1.0);
                
                RCP<OperatorT<MultiVector<SC,LO,GO,NO> > > OpK = rcp(new XpetraOp<SC, LO, GO, NO>(K1));
                RCP<OperatorT<MultiVector<SC,LO,GO,NO> > > OpP = rcp(new XpetraOp<SC, LO, GO, NO>(OneLevelPrec));
                
                RCP<LinearProblem<SC,MultiVector<SC,LO,GO,NO>,OperatorT<MultiVector<SC,LO,GO,NO> > > > belosLinearProblem(new LinearProblem<SC,MultiVector<SC,LO,GO,NO>,OperatorT<MultiVector<SC,LO,GO,NO> > >(OpK,xSolution,xRightHandSide));
                SolverFactory<SC,MultiVector<SC,LO,GO,NO>,OperatorT<MultiVector<SC,LO,GO,NO> > > belosFactory;
                RCP<ParameterList> solverParameterList = sublist(parameterList,"Solver");
                RCP<SolverManager<SC,MultiVector<SC,LO,GO,NO>,OperatorT<MultiVector<SC,LO,GO,NO> > > > belosSoverManager = belosFactory.create(solverParameterList->get("Solver","GMRES"),sublist(solverParameterList,solverParameterList->get("Solver","GMRES")));
                belosSoverManager->setProblem(belosLinearProblem);
                
                if (!solverParameterList->get("PreconditionerPosition","left").compare("left")) {
                    belosLinearProblem->setLeftPrec(OpP);
                } else if (!solverParameterList->get("PreconditionerPosition","left").compare("right")) {
                    belosLinearProblem->setRightPrec(OpP);
                } else {
                    FROSCH_ASSERT(0!=0,"PreconditionerPosition unknown...");
                }
                
                belosLinearProblem->setProblem(xSolution,xRightHandSide);
                belosSoverManager->solve();
                belosSoverManager->getNumIters();

            }
            else{
                
                RCP<TwoLevelBlockPreconditioner<SC,LO,GO,NO> > TwoLevelPrec(new TwoLevelBlockPreconditioner<SC,LO,GO,NO>(KFull,sublist(parameterList,"TwoLevelPreconditioner")));
                //
                if (Comm->MyPID()==0) cout << "INITIALIZE...";
                SchwarzOperator<SC,LO,GO,NO>::GOVecPtr blockMaxGIDVector(2);
                blockMaxGIDVector[0] = K1->getRowMap()->getMaxAllGlobalIndex();
                blockMaxGIDVector[1] = mapFull->getMaxAllGlobalIndex();
                
                SchwarzOperator<SC,LO,GO,NO>::UNVecPtr dofsPerNodeVec(2);
                Teuchos::ArrayRCP<DofOrdering> ordVec(2);
                dofsPerNodeVec[0] = 1;
                dofsPerNodeVec[1] = 1;
                ordVec[0] = NodeWise;
                ordVec[1] = NodeWise;
                TwoLevelPrec->initialize(Dimension,dofsPerNodeVec,ordVec,blockMaxGIDVector,Overlap);
                
                if (Comm->MyPID()==0) cout << "COMPUTE...";
                TwoLevelPrec->compute();
                if (Comm->MyPID()==0) cout << "done" << endl << "SOLVING EQUATION SYSTEM...";
                
                RCP<MultiVector<SC,LO,GO,NO> > xSolution = MultiVectorFactory<SC,LO,GO,NO>::Build(mapFull,1);
                RCP<MultiVector<SC,LO,GO,NO> > xRightHandSide = MultiVectorFactory<SC,LO,GO,NO>::Build(mapFull,1);
                
                xSolution->putScalar(0.0);
                xRightHandSide->putScalar(1.0);
                
                RCP<OperatorT<MultiVector<SC,LO,GO,NO> > > OpK = rcp(new XpetraOp<SC, LO, GO, NO>(KFull));
                RCP<OperatorT<MultiVector<SC,LO,GO,NO> > > OpP = rcp(new XpetraOp<SC, LO, GO, NO>(TwoLevelPrec));
                
                RCP<LinearProblem<SC,MultiVector<SC,LO,GO,NO>,OperatorT<MultiVector<SC,LO,GO,NO> > > > belosLinearProblem(new LinearProblem<SC,MultiVector<SC,LO,GO,NO>,OperatorT<MultiVector<SC,LO,GO,NO> > >(OpK,xSolution,xRightHandSide));
                SolverFactory<SC,MultiVector<SC,LO,GO,NO>,OperatorT<MultiVector<SC,LO,GO,NO> > > belosFactory;
                RCP<ParameterList> solverParameterList = sublist(parameterList,"Solver");
                RCP<SolverManager<SC,MultiVector<SC,LO,GO,NO>,OperatorT<MultiVector<SC,LO,GO,NO> > > > belosSoverManager = belosFactory.create(solverParameterList->get("Solver","GMRES"),sublist(solverParameterList,solverParameterList->get("Solver","GMRES")));
                belosSoverManager->setProblem(belosLinearProblem);
                
                if (!solverParameterList->get("PreconditionerPosition","left").compare("left")) {
                    belosLinearProblem->setLeftPrec(OpP);
                } else if (!solverParameterList->get("PreconditionerPosition","left").compare("right")) {
                    belosLinearProblem->setRightPrec(OpP);
                } else {
                    FROSCH_ASSERT(0!=0,"PreconditionerPosition unknown...");
                }
                
                belosLinearProblem->setProblem(xSolution,xRightHandSide);
                belosSoverManager->solve();
                belosSoverManager->getNumIters();

            }
            

//			RCP<TwoLevelBlockPreconditioner<SC,LO,GO,NO> > TwoLevelPrec(new TwoLevelBlockPreconditioner<SC,LO,GO,NO>(KFull,sublist(parameterList,"TwoLevelPreconditioner")));
////
//			if (Comm->MyPID()==0) cout << "INITIALIZE...";
//            SchwarzOperator<SC,LO,GO,NO>::GOVecPtr blockMaxGIDVector(2);
//            blockMaxGIDVector[0] = K1->getRowMap()->getMaxAllGlobalIndex();
//            blockMaxGIDVector[1] = mapFull->getMaxAllGlobalIndex();
//            
//            SchwarzOperator<SC,LO,GO,NO>::UNVecPtr dofsPerNodeVec(2);
//            Teuchos::ArrayRCP<DofOrdering> ordVec(2);
//            dofsPerNodeVec[0] = 1;
//            dofsPerNodeVec[1] = 1;
//            ordVec[0] = NodeWise;
//            ordVec[1] = NodeWise;
//			TwoLevelPrec->initialize(Dimension,Overlap,dofsPerNodeVec,ordVec,blockMaxGIDVector);
//            
//			if (Comm->MyPID()==0) cout << "COMPUTE...";
//			TwoLevelPrec->compute();
//			if (Comm->MyPID()==0) cout << "done" << endl << "SOLVING EQUATION SYSTEM...";
//
//			RCP<MultiVector<SC,LO,GO,NO> > xSolution = MultiVectorFactory<SC,LO,GO,NO>::Build(mapFull,1);
//			RCP<MultiVector<SC,LO,GO,NO> > xRightHandSide = MultiVectorFactory<SC,LO,GO,NO>::Build(mapFull,1);
//
//			xSolution->putScalar(0.0);
//			xRightHandSide->putScalar(1.0);
//
//			RCP<OperatorT<MultiVector<SC,LO,GO,NO> > > OpK = rcp(new XpetraOp<SC, LO, GO, NO>(KFull));
////			RCP<OperatorT<MultiVector<SC,LO,GO,NO> > > OpP = rcp(new XpetraOp<SC, LO, GO, NO>(TwoLevelPrec));
//            RCP<OperatorT<MultiVector<SC,LO,GO,NO> > > OpP = rcp(new XpetraOp<SC, LO, GO, NO>(OneLevelPrec));
//
//			RCP<LinearProblem<SC,MultiVector<SC,LO,GO,NO>,OperatorT<MultiVector<SC,LO,GO,NO> > > > belosLinearProblem(new LinearProblem<SC,MultiVector<SC,LO,GO,NO>,OperatorT<MultiVector<SC,LO,GO,NO> > >(OpK,xSolution,xRightHandSide));
//			SolverFactory<SC,MultiVector<SC,LO,GO,NO>,OperatorT<MultiVector<SC,LO,GO,NO> > > belosFactory;
//			RCP<ParameterList> solverParameterList = sublist(parameterList,"Solver");
//			RCP<SolverManager<SC,MultiVector<SC,LO,GO,NO>,OperatorT<MultiVector<SC,LO,GO,NO> > > > belosSoverManager = belosFactory.create(solverParameterList->get("Solver","GMRES"),sublist(solverParameterList,solverParameterList->get("Solver","GMRES")));
//			belosSoverManager->setProblem(belosLinearProblem);
//
//			if (!solverParameterList->get("PreconditionerPosition","left").compare("left")) {
//				belosLinearProblem->setLeftPrec(OpP);
//			} else if (!solverParameterList->get("PreconditionerPosition","left").compare("right")) {
//				belosLinearProblem->setRightPrec(OpP);
//			} else {
//				FROSCH_ASSERT(0!=0,"PreconditionerPosition unknown...");
//			}
//
//			belosLinearProblem->setProblem(xSolution,xRightHandSide);
//			belosSoverManager->solve();
//			belosSoverManager->getNumIters();

			if (Comm->MyPID()==0) cout << "done" << endl;
		}

		MPI_Comm_free(&COMM);

	}

	MPI_Finalize();

	return(EXIT_SUCCESS);
}
