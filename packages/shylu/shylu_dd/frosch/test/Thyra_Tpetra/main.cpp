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

#include <FROSch_GDSWPreconditioner_def.hpp>
#include <FROSch_RGDSWPreconditioner_def.hpp>

// Thyra includes
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_SolveSupportTypes.hpp>
#include <Thyra_BelosLinearOpWithSolveFactory.hpp>
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>
#include <Thyra_TpetraLinearOp.hpp>
#include <Thyra_TpetraMultiVector.hpp>
#include <Thyra_TpetraVector.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_VectorStdOps.hpp>
#include <Thyra_EpetraLinearOp.hpp>
#include <Thyra_VectorSpaceBase_def.hpp>
#include <Thyra_VectorSpaceBase_decl.hpp>
#include <Thyra_BelosLinearOpWithSolve_def.hpp>
#include <Thyra_BelosLinearOpWithSolveFactory_def.hpp>

// Stratimikos includes
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Stratimikos_FROSchXpetra.hpp>

// Xpetra include
#include <Xpetra_Parameters.hpp>

// FROSCH thyra includes
#include "Thyra_FROSchLinearOp_def.hpp"
#include "Thyra_FROSchXpetraFactory_def.hpp"
#include <MueLu_Utilities_def.hpp>
#include <MueLu_Utilities_decl.hpp>



typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef int GO;
typedef KokkosClassic::DefaultNode::DefaultNodeType EpetraNode;
typedef EpetraNode NO;
typedef KokkosClassic::DefaultNode::DefaultNodeType TpetraNode;
typedef TpetraNode TO;

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
        Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();


		int M = 4;
		My_CLP.setOption("M",&M,"H / h.");
		int Dimension = 2;
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
				parameterList = getParametersFromXmlFile("xpetra_ParameterList.xml");
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
			//        RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(std::cout)); K->describe(*fancy,Teuchos::VERB_EXTREME);
    
            //RCP<Matrix<SC,LO,GO,TO> > TK = rcp_dynamic_cast<Matrix<SC,LO,GO,TO> >(K);
            RCP<const Xpetra::CrsMatrixWrap<SC,LO,GO,NO> > crsOp = rcp_dynamic_cast<const Xpetra::CrsMatrixWrap<SC,LO,GO,NO> >(K);
            if (crsOp == Teuchos::null){
                if (Comm->MyPID()==0) cout << "CRS Cast not working\n";
                
                }
            const RCP<const Xpetra::TpetraCrsMatrix<SC,LO,GO,NO> > &tmp_ECrsMtx = rcp_dynamic_cast<const Xpetra::TpetraCrsMatrix<SC,LO,GO,NO> >(crsOp->getCrsMatrix());
            if (tmp_ECrsMtx == Teuchos::null){
                if (Comm->MyPID()==0) cout << "CRS Tpetra Cast not working\n";

            }
            RCP<Tpetra::CrsMatrix<SC,LO,GO,NO> > TMat = tmp_ECrsMtx->getTpetra_CrsMatrixNonConst();

            //TEUCHOS_TEST_FOR_EXCEPT(Teuchos::is_null(TK));

			RCP<MultiVector<SC,LO,GO,NO> > xSolution = MultiVectorFactory<SC,LO,GO,NO>::Build(UniqueMap,1);
			RCP<MultiVector<SC,LO,GO,NO> > xRightHandSide = MultiVectorFactory<SC,LO,GO,NO>::Build(UniqueMap,1);

			xSolution->putScalar(0.0);
            xRightHandSide->putScalar(1.0);
            
            //Teuchos::RCP<Xpetra::CrsMatrix<SC, LO, GO, NO > > exA = Teuchos::rcp(new Xpetra::TpetraCrsMatrix<int,TO>(K));
            //Teuchos::RCP<CrsMatrixWrap<SC,LO,GO> > exAWrap = Teuchos::rcp(new CrsMatrixWrap<SC,LO,GO> (exA));
            
            RCP<const Thyra::LinearOpBase<SC> > K_thyra = Xpetra::ThyraUtils<SC,LO,GO,NO>::toThyra(crsOp->getCrsMatrix());
            RCP<Thyra::MultiVectorBase<SC> >thyraX =
            Teuchos::rcp_const_cast<Thyra::MultiVectorBase<SC> >(Xpetra::ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(xSolution));
            RCP<const Thyra::MultiVectorBase<SC> >thyraB = Xpetra::ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(xRightHandSide);
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "----------------done-----------\n";
            
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "----------------Stratimikos LinearSolverBuilder-----------\n";
            
            //-----------Set Coordinates and RepMap in ParameterList--------------------------
            // RCP<ParameterList> plList =  sublist(parameterList,"Preconditioner Types");
            // sublist(plList,"TwoLevelPreconditioner")->set("Coordinates",Coord);
            // sublist(plList,"TwoLevelPreconditioner")->set("RepeatedMap",RepMapX);
            
            Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
            Stratimikos::enableFROSch<LO,GO,NO>(linearSolverBuilder);
            
            linearSolverBuilder.setParameterList(parameterList);
            
            
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "----------------done-----------\n";
            
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "----------------Thyra PrepForSolve-----------\n";
            
            
            
            RCP<Thyra::LinearOpWithSolveFactoryBase<SC> > lowsFactory =
            linearSolverBuilder.createLinearSolveStrategy("");
            
            lowsFactory->setOStream(out);
            lowsFactory->setVerbLevel(Teuchos::VERB_HIGH);
            
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "----------------done-----------\n";
            if (Comm->MyPID()==0) cout << "----------------Thyra LinearOpWithSolve-----------\n";
            
            RCP<Thyra::LinearOpWithSolveBase<SC> > lows =
            Thyra::linearOpWithSolve(*lowsFactory, K_thyra);
            
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "----------------done-----------\n";
            
            if (Comm->MyPID()==0) cout << "----------------Solve-----------\n";
            Thyra::SolveStatus<double> status =
            Thyra::solve<double>(*lows, Thyra::NOTRANS, *thyraB, thyraX.ptr());
            
            
            Comm->Barrier();
            if (Comm->MyPID()==0) cout << "----------------done-----------\n";
            

            
		}

		MPI_Comm_free(&COMM);

	}

	MPI_Finalize();

	return(EXIT_SUCCESS);
}
