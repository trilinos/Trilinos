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

#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_MultiVectorIn.h"

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

#include <FROSch_TwoLevelPreconditioner_def.hpp>

//#include "Tools/FROSch_Tools_def.hpp"

typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef int GO;
typedef KokkosClassic::DefaultNode::DefaultNodeType EpetraNode;
typedef EpetraNode NO;

using namespace std;
using namespace Teuchos;
using namespace Xpetra;
using namespace FROSch;
using namespace Belos;
using namespace EpetraExt;

int main(int argc, char *argv[])
{
	MPI_Init(&argc,&argv);

	{

		RCP<Epetra_MpiComm> Comm(new Epetra_MpiComm(MPI_COMM_WORLD));
		RCP<const Teuchos::Comm<int> > TeuchosComm = rcp(new MpiComm<int> (MPI_COMM_WORLD));

		RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout));

		{
			///////////////////
			// ParameterList //
			///////////////////
			RCP<ParameterList> parameterList = getParametersFromXmlFile("Parameters.xml");
			if (Comm->MyPID()==0) {
				cout << "--------------------------------------------------------------------------------\nPARAMETERS:" << endl;
				parameterList->print(cout);
				cout << "--------------------------------------------------------------------------------\n\n";
			}

			//////////////////
			// Repeated Map //
			//////////////////
			Epetra_Map *repeatedMap;
			MatrixMarketFileToMap("repeatedMap.dat",*Comm,repeatedMap);
			RCP<Map<LO,GO,NO> > RepeatedMap = ConvertToXpetra<LO,GO,NO>(UseTpetra,*repeatedMap,TeuchosComm);

			////////////////
			// Unique Map //
			////////////////
			Teuchos::RCP<Xpetra::Map<LO,GO,NO> > repeatedNodesMap;
			Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::Map<LO,GO,NO> > > dofsMaps;
			BuildDofMaps(RepeatedMap,2,NodeWise,repeatedNodesMap,dofsMaps);
			RCP<Map<LO,GO,NO> > uniqueNodesMap = BuildUniqueMap<LO,GO,NO>(repeatedNodesMap);
			RCP<Map<LO,GO,NO> > UniqueMap = BuildMapFromNodeMap(uniqueNodesMap,2,NodeWise);

			//////////////////////
			// Stiffness Matrix //
			//////////////////////
			Epetra_CrsMatrix *KEpetra;
			MatrixMarketFileToCrsMatrix("system_matrix.dat",*Comm,KEpetra);
			RCP<Matrix<SC,LO,GO,NO> > Ktmp = ConvertToXpetra<SC,LO,GO,NO>(UseTpetra,*KEpetra,TeuchosComm);
			//        Ktmp->describe(*fancy,VERB_EXTREME);

			RCP<Matrix<SC,LO,GO,NO> > K = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(UniqueMap,2*Ktmp->getGlobalMaxNumRowEntries());
			RCP<Import<LO,GO,NO> > scatter = Xpetra::ImportFactory<LO,GO,NO>::Build(Ktmp->getRowMap(),UniqueMap);
			K->doImport(*Ktmp,*scatter,Xpetra::ADD);
			K->fillComplete();
			//    K->describe(*fancy,VERB_EXTREME);

			///////////////
			// Node List //
			///////////////
			Epetra_Map *nodeMapEpetra;
			EpetraExt::MatrixMarketFileToMap("nodeMap.dat",*Comm,nodeMapEpetra);
			Epetra_MultiVector *nodeList;
			EpetraExt::MatrixMarketFileToMultiVector("nodeList.dat",*nodeMapEpetra,nodeList);
			RCP<MultiVector<SC,LO,GO,NO> > NodeList = ConvertToXpetra<SC,LO,GO,NO>(UseTpetra,*nodeList,TeuchosComm);
			//        NodeList->describe(*fancy,VERB_EXTREME);


			if (Comm->MyPID()==0) cout << "CONSTRUCTING PRECONDITIONER...";
			RCP<TwoLevelPreconditioner<SC,LO,GO,NO> > TwoLevelPrec(new TwoLevelPreconditioner<SC,LO,GO,NO>(K,sublist(parameterList,"TwoLevelPreconditioner")));

			if (Comm->MyPID()==0) cout << "INITIALIZE...";
			TwoLevelPrec->initialize(2,1,RepeatedMap,2,NodeWise,NodeList);
			if (Comm->MyPID()==0) cout << "COMPUTE...";
			TwoLevelPrec->compute();
			if (Comm->MyPID()==0) cout << "done" << endl << "SOLVING EQUATION SYSTEM...";
			RCP<MultiVector<SC,LO,GO,NO> > xSolution = MultiVectorFactory<SC,LO,GO,NO>::Build(UniqueMap,1);
			RCP<MultiVector<SC,LO,GO,NO> > xRightHandSide = MultiVectorFactory<SC,LO,GO,NO>::Build(UniqueMap,1);

			xSolution->putScalar(0.0);
			xRightHandSide->putScalar(1.0);

			RCP<OperatorT<MultiVector<SC,LO,GO,NO> > > OpK = rcp(new XpetraOp<SC, LO, GO, NO>(K));
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

			if (Comm->MyPID()==0) cout << "done" << endl;

		}

	}

	MPI_Finalize();

	return(EXIT_SUCCESS);
}
