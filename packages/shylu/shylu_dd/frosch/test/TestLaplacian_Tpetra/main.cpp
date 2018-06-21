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
#include <BelosSolverFactory_Generic.hpp>
//#include <BelosPseudoBlockGmresSolMgr.hpp>

#include <FROSch_GDSWPreconditioner_def.hpp>
#include <FROSch_RGDSWPreconditioner_def.hpp>

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

int main(int argc, char *argv[])
{
    MPI_Init(&argc,&argv);
    Epetra_MpiComm CommWorld(MPI_COMM_WORLD);
    
    CommandLineProcessor My_CLP;
    
    int M = 4;
    My_CLP.setOption("M",&M,"H / h.");
    int Dimension = 3;
    My_CLP.setOption("DIM",&Dimension,"Dimension.");
    int Overlap = 1;
    My_CLP.setOption("OL",&Overlap,"Overlap.");
    bool Reduced = false;
    My_CLP.setOption("RGDSW","GDSW",&Reduced,"Using the reduced coarse space.");
    
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
        
        RCP<ParameterList> parameterList;
        if (!Reduced) {
            parameterList = getParametersFromXmlFile("ParametersGDSW.xml");
        } else {
            parameterList = getParametersFromXmlFile("ParametersRGDSW.xml");
        }
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
        RCP<Epetra_CrsMatrix> KEpetra;
        
        if (Dimension==2) {
            UniqueMapEpetra.reset(Galeri::CreateMap("Cartesian2D", *Comm, GalerList));
            KEpetra.reset(Galeri::CreateCrsMatrix("Laplace2D", UniqueMapEpetra.get(), GalerList));
        } else if (Dimension==3) {
            UniqueMapEpetra.reset(Galeri::CreateMap("Cartesian3D", *Comm, GalerList));
            KEpetra.reset(Galeri::CreateCrsMatrix("Laplace3D", UniqueMapEpetra.get(), GalerList));
        }
        
        ArrayView<GO> uniqueMapArrayView(UniqueMapEpetra->MyGlobalElements(),UniqueMapEpetra->NumMyElements());
        RCP<Map<LO,GO,NO> > UniqueMap = MapFactory<LO,GO,NO>::Build(UseTpetra,-1,uniqueMapArrayView,0,TeuchosComm);
        RCP<Matrix<SC,LO,GO,NO> > K = MatrixFactory<SC,LO,GO,NO>::Build(UniqueMap,KEpetra->MaxNumEntries());
        for (LO i=0; i<UniqueMapEpetra->NumMyElements(); i++) {
            LO numEntries;
            GO* indices;
            SC* values;
            KEpetra->ExtractMyRowView(i,numEntries,values,indices);
            
            Array<GO> indicesArray(numEntries);
            ArrayView<SC> valuesArrayView(values,numEntries);
            for (LO j=0; j<numEntries; j++) {
                indicesArray[j] = KEpetra->ColMap().GID(indices[j]);
            }
            K->insertGlobalValues(UniqueMapEpetra->GID(i),indicesArray(),valuesArrayView);
        }
        K->fillComplete();
        
        if (Comm->MyPID()==0) cout << "done" << endl << "CONSTRUCTING PRECONDITIONER...";
        RCP<SchwarzPreconditioner<SC,LO,GO,NO> > Preconditioner;
        if (!Reduced) {
            RCP<GDSWPreconditioner<SC,LO,GO,NO> > TmpPrec(new GDSWPreconditioner<SC,LO,GO,NO>(K,sublist(parameterList,"GDSWPreconditioner")));
            if (Comm->MyPID()==0) cout << "INITIALIZE...";
            TmpPrec->initialize(Dimension,1);
            Preconditioner = TmpPrec;
        } else {
            RCP<RGDSWPreconditioner<SC,LO,GO,NO> > TmpPrec(new RGDSWPreconditioner<SC,LO,GO,NO>(K,sublist(parameterList,"RGDSWPreconditioner")));
            if (Comm->MyPID()==0) cout << "INITIALIZE...";
            TmpPrec->initialize(Dimension,1);
            Preconditioner = TmpPrec;
        }
        
        if (Comm->MyPID()==0) cout << "COMPUTE...";
        Preconditioner->compute();
        if (Comm->MyPID()==0) cout << "done" << endl << "SOLVING EQUATION SYSTEM...";
        
        RCP<MultiVector<SC,LO,GO,NO> > xSolution = MultiVectorFactory<SC,LO,GO,NO>::Build(UniqueMap,1);
        RCP<MultiVector<SC,LO,GO,NO> > xRightHandSide = MultiVectorFactory<SC,LO,GO,NO>::Build(UniqueMap,1);
        
        xSolution->putScalar(0.0);
        xRightHandSide->putScalar(1.0);
        
        RCP<OperatorT<MultiVector<SC,LO,GO,NO> > > OpK = rcp(new XpetraOp<SC, LO, GO, NO>(K));
        RCP<OperatorT<MultiVector<SC,LO,GO,NO> > > OpP = rcp(new XpetraOp<SC, LO, GO, NO>(Preconditioner));
        
        RCP<LinearProblem<SC,MultiVector<SC,LO,GO,NO>,OperatorT<MultiVector<SC,LO,GO,NO> > > > belosLinearProblem(new LinearProblem<SC,MultiVector<SC,LO,GO,NO>,OperatorT<MultiVector<SC,LO,GO,NO> > >(OpK,xSolution,xRightHandSide));
        GenericSolverFactory<SC,MultiVector<SC,LO,GO,NO>,OperatorT<MultiVector<SC,LO,GO,NO> > > belosFactory;

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
    
    MPI_Finalize();
    
    return(EXIT_SUCCESS);
}
