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
#define FROSCH_OFFSET_MAPS

#include <mpi.h>
#include <Epetra_MpiComm.h>

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#if defined(HAVE_SHYLU_DDFROSCH_EPETRAEXT) && defined(HAVE_EPETRAEXT_HDF5)
#include "EpetraExt_BlockMapIn.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_MultiVectorIn.h"
#include "EpetraExt_HDF5.h"
#endif
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
#if defined(HAVE_SHYLU_DDFROSCH_EPETRAEXT) && defined(HAVE_EPETRAEXT_HDF5)
            unsigned Dimension = 2;
            RCP<HDF5> hDF5IO(new HDF5(*Comm));
            hDF5IO->Open("stokes.h5");
            string groupNameMatrix = "Matrix";
            string groupNameRepeatedMapVelo =  "RepeatedMapVelocity";
            string groupNameRepeatedMapPress =  "RepeatedMapPressure";
            string groupNameRHS = "RHS";
            
            
            ///////////////////
            // Repeated Maps //
            ///////////////////
            Epetra_Map *repeatedMapEpetraVelo;
            hDF5IO->Read(groupNameRepeatedMapVelo,repeatedMapEpetraVelo);
            RCP<Map<LO,GO,NO> > repeatedMapVelo = ConvertToXpetra<LO,GO,NO>(UseTpetra,*repeatedMapEpetraVelo,TeuchosComm);

            Epetra_Map *repeatedMapEpetraPress;
            hDF5IO->Read(groupNameRepeatedMapPress,repeatedMapEpetraPress);
            GO offsetVelocityMap = repeatedMapVelo->getMaxAllGlobalIndex()+1;
            Array<GO> elementList(repeatedMapEpetraPress->NumMyElements());
            for (unsigned i=0; i<elementList.size(); i++) {
                elementList[i] = repeatedMapEpetraPress->GID(i) + offsetVelocityMap;
            }
            RCP<Map<LO,GO,NO> > repeatedMapPress = Xpetra::MapFactory<LO,GO,NO>::Build(UseTpetra,-1,elementList,0,TeuchosComm);
            
            Teuchos::ArrayRCP<Teuchos::RCP<Xpetra::Map<LO,GO,NO> > > repeatedMapsVector(2);

            repeatedMapsVector[0] = repeatedMapVelo;
            repeatedMapsVector[1] = repeatedMapPress;
            
            RCP<Map<LO,GO,NO> > repeatedMap = MergeMaps(repeatedMapsVector);

            ////////////////
            // Unique Map //
            ////////////////
            RCP<Map<LO,GO,NO> > uniqueMap = BuildUniqueMap<LO,GO,NO>(repeatedMap);
            
            
            /////////
            // RHS //
            /////////
            Epetra_MultiVector *rhsEpetra;
            hDF5IO->Read(groupNameRHS,rhsEpetra);

            RCP<MultiVector<SC,LO,GO,NO> > rhsTmp = ConvertToXpetra<SC,LO,GO,NO>(UseTpetra,*rhsEpetra,TeuchosComm);
            RCP<MultiVector<SC,LO,GO,NO> > rhs = Xpetra::MultiVectorFactory<SC,LO,GO,NO>::Build(uniqueMap,1);
            RCP<Import<LO,GO,NO> > scatter = Xpetra::ImportFactory<LO,GO,NO>::Build(rhsTmp->getMap(),uniqueMap);
            rhs->doImport(*rhsTmp,*scatter,Xpetra::ADD);
            
            ////////////
            // Matrix //
            ////////////
            Epetra_CrsMatrix *matrixEpetra;
            hDF5IO->Read(groupNameMatrix,matrixEpetra);
            RCP<Matrix<SC,LO,GO,NO> > matrixTmp = ConvertToXpetra<SC,LO,GO,NO>(UseTpetra,*matrixEpetra,TeuchosComm);
            RCP<Matrix<SC,LO,GO,NO> > matrix = Xpetra::MatrixFactory<SC,LO,GO,NO>::Build(uniqueMap,matrixTmp->getGlobalMaxNumRowEntries());

            matrix->doImport(*matrixTmp,*scatter,Xpetra::ADD);
            matrix->fillComplete();
            
			if (Comm->MyPID()==0) cout << "CONSTRUCTING PRECONDITIONER...";
            
        
            RCP<TwoLevelBlockPreconditioner<SC,LO,GO,NO> > TwoLevelPrec(new TwoLevelBlockPreconditioner<SC,LO,GO,NO>(matrix,sublist(parameterList,"TwoLevelPreconditioner")));

            if (Comm->MyPID()==0) cout << "INITIALIZE...";
            SchwarzOperator<SC,LO,GO,NO>::GOVecPtr blockMaxGID(2);
            blockMaxGID[0] = repeatedMapVelo->getMaxAllGlobalIndex(); //CH: 30.08.18 fix offet, not needed anymore? Only needed for setup without repeated Maps
            blockMaxGID[1] = repeatedMapPress->getMaxAllGlobalIndex();
            
            Teuchos::ArrayRCP<unsigned> dofsPerNodeVec(2);
            dofsPerNodeVec[0] = 2;
            dofsPerNodeVec[1] = 1;
            Teuchos::ArrayRCP<DofOrdering> ordVec(2);
            ordVec[0] = NodeWise;
            ordVec[1] = NodeWise;
            
            int Overlap = parameterList->sublist("TwoLevelPreconditioner").get("Overlap",1);
            TwoLevelPrec->initialize(Dimension,dofsPerNodeVec,ordVec,blockMaxGID,Overlap,repeatedMapsVector);
        
            if (Comm->MyPID()==0) cout << "COMPUTE...";
            TwoLevelPrec->compute();
            
            if (Comm->MyPID()==0) cout << "done" << endl << "SOLVING EQUATION SYSTEM...";
            RCP<MultiVector<SC,LO,GO,NO> > xSolution = MultiVectorFactory<SC,LO,GO,NO>::Build(uniqueMap,1);

            xSolution->putScalar(0.0);

            RCP<OperatorT<MultiVector<SC,LO,GO,NO> > > OpK = rcp(new XpetraOp<SC, LO, GO, NO>(matrix));
            RCP<OperatorT<MultiVector<SC,LO,GO,NO> > > OpP = rcp(new XpetraOp<SC, LO, GO, NO>(TwoLevelPrec));

            RCP<LinearProblem<SC,MultiVector<SC,LO,GO,NO>,OperatorT<MultiVector<SC,LO,GO,NO> > > > belosLinearProblem(new LinearProblem<SC,MultiVector<SC,LO,GO,NO>,OperatorT<MultiVector<SC,LO,GO,NO> > >(OpK,xSolution,rhs));
            SolverFactory<SC,MultiVector<SC,LO,GO,NO>,OperatorT<MultiVector<SC,LO,GO,NO> > > belosFactory;
            RCP<ParameterList> solverParameterList = sublist(parameterList,"Solver");
            RCP<SolverManager<SC,MultiVector<SC,LO,GO,NO>,OperatorT<MultiVector<SC,LO,GO,NO> > > > belosSoverManager = belosFactory.create(solverParameterList->get("Solver","GMRES"),sublist(solverParameterList,solverParameterList->get("Solver","GMRES")));
            belosSoverManager->setProblem(belosLinearProblem);

            if (!solverParameterList->get("PreconditionerPosition","left").compare("left")) {
                belosLinearProblem->setLeftPrec(OpP);
            } else if (!solverParameterList->get("PreconditionerPosition","left").compare("right")) {
                belosLinearProblem->setRightPrec(OpP);
            } else {
                FROSCH_ASSERT(false,"PreconditionerPosition unknown...");
            }
            if (!parameterList->sublist("TwoLevelPreconditioner").get("Level Combination","Additive").compare("Multiplicative")) {
                TwoLevelPrec->preApplyCoarse(rhs,xSolution);
            }
            belosLinearProblem->setProblem(xSolution,rhs);
            belosSoverManager->solve();
            belosSoverManager->getNumIters();
            
			if (Comm->MyPID()==0) cout << "done" << endl;

#endif
		}

	}

	MPI_Finalize();

	return(EXIT_SUCCESS);
}
