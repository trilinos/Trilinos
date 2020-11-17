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

#include <ShyLU_DDFROSch_config.h>

#include <mpi.h>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>
#include <Teuchos_StackedTimer.hpp>

#include <Epetra_MpiComm.h>
#if defined(HAVE_SHYLU_DDFROSCH_EPETRAEXT) && defined(HAVE_SHYLU_DDFROSCH_HDF5)
#include "EpetraExt_HDF5.h"
#endif

// Galeri::Xpetra
#include "Galeri_XpetraProblemFactory.hpp"
#include "Galeri_XpetraMatrixTypes.hpp"
#include "Galeri_XpetraParameters.hpp"
#include "Galeri_XpetraUtils.hpp"
#include "Galeri_XpetraMaps.hpp"

// Thyra includes
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_SolveSupportTypes.hpp>
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>
#include <Thyra_TpetraLinearOp.hpp>
#include <Thyra_TpetraMultiVector.hpp>
#include <Thyra_TpetraVector.hpp>
#include <Thyra_TpetraThyraWrappers.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_VectorStdOps.hpp>
#ifdef HAVE_SHYLU_DDFROSCH_EPETRA
#include <Thyra_EpetraLinearOp.hpp>
#endif
#include <Thyra_VectorSpaceBase_def.hpp>
#include <Thyra_VectorSpaceBase_decl.hpp>

// Stratimikos includes
//#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Stratimikos_FROSchXpetra.hpp>

#include <Tpetra_Core.hpp>

// Xpetra include
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_DefaultPlatform.hpp>
#ifdef HAVE_XPETRA_EPETRA
#include <Xpetra_EpetraCrsMatrix.hpp>
#endif
#include <Xpetra_Parameters.hpp>

// FROSCH thyra includes
#include "Thyra_FROSchLinearOp_def.hpp"
#include "Thyra_FROSchFactory_def.hpp"
#include <FROSch_Tools_def.hpp>


using UN    = unsigned;
using SC    = double;
using LO    = int;
using GO    = FROSch::DefaultGlobalOrdinal;
using NO    = KokkosClassic::DefaultNode::DefaultNodeType;

using namespace std;
using namespace EpetraExt;
using namespace FROSch;
using namespace Teuchos;
using namespace Xpetra;
using namespace Thyra;

int main(int argc, char *argv[])
{
    oblackholestream blackhole;
    GlobalMPISession mpiSession(&argc,&argv,&blackhole);

    RCP<const Comm<int> > CommWorld = DefaultPlatform::getDefaultPlatform().getComm();

    CommandLineProcessor My_CLP;

    RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();

    // Overlap
    int Overlap = 0;
    My_CLP.setOption("O",&Overlap,"Overlap.");

    // XML-file specifying the solver configuration
    string xmlFile = "ParameterList_TwoLevelBlockPreconditioner_GDSW.xml";
    My_CLP.setOption("PLIST",&xmlFile,"File name of the parameter list.");

    // Use Epetra as linear algebra stack
    bool useepetra = false;
    My_CLP.setOption("USEEPETRA","USETPETRA",&useepetra,"Use Epetra infrastructure for the linear algebra.");

    // File name of HDF5 file to be read
    string hdf5File = "stokes.h5";
    My_CLP.setOption("HDF5FILE", &hdf5File, "File name of HDF5 file to be read.");

    My_CLP.recogniseAllOptions(true);
    My_CLP.throwExceptions(false);
    CommandLineProcessor::EParseCommandLineReturn parseReturn = My_CLP.parse(argc,argv);
    if (parseReturn == CommandLineProcessor::PARSE_HELP_PRINTED) {
        return(EXIT_SUCCESS);
    }

    CommWorld->barrier();
    RCP<StackedTimer> stackedTimer = rcp(new StackedTimer("Thyra Stokes Test"));
    TimeMonitor::setStackedTimer(stackedTimer);

    int color=1;
    if (CommWorld->getRank()<4) {
        color=0;
    }

    UnderlyingLib xpetraLib = UseTpetra;
    if (useepetra) {
        xpetraLib = UseEpetra;
    } else {
        xpetraLib = UseTpetra;
    }

    RCP<const Comm<int> > Comm = CommWorld->split(color,CommWorld->getRank());

    MPI_Comm COMM;
    MPI_Comm_split(MPI_COMM_WORLD,color,CommWorld->getRank(),&COMM);
    RCP<Epetra_MpiComm> EpetraComm(new Epetra_MpiComm(COMM));

#ifdef HAVE_EPETRAEXT_HDF5
    if (color==0) {

        RCP<ParameterList> parameterList = getParametersFromXmlFile(xmlFile);

        Comm->barrier(); if (Comm->getRank()==0) cout << "##############################\n# Import Monolythic System #\n##############################\n" << endl;

        unsigned Dimension = 2;
        RCP<HDF5> hDF5IO(new HDF5(*EpetraComm));
        hDF5IO->Open(hdf5File);

        ///////////////////
        // Repeated Maps //
        ///////////////////
        string groupNameRepeatedMapVelo =  "RepeatedMapVelocity";
        Epetra_Map *repeatedMapEpetraVelo;
        hDF5IO->Read(groupNameRepeatedMapVelo,repeatedMapEpetraVelo);

        RCP<Map<LO,GO,NO> > repeatedMapVelo = ConvertToXpetra<SC,LO,GO,NO>::ConvertMap(xpetraLib,*repeatedMapEpetraVelo,Comm);

        string groupNameRepeatedMapPress =  "RepeatedMapPressure";
        Epetra_Map *repeatedMapEpetraPress;
        hDF5IO->Read(groupNameRepeatedMapPress,repeatedMapEpetraPress);
        RCP<Map<LO,GO,NO> > repeatedMapPress = ConvertToXpetra<SC,LO,GO,NO>::ConvertMap(xpetraLib,*repeatedMapEpetraPress,Comm);

//        GO offsetVelocityMap = repeatedMapVelo->getMaxAllGlobalIndex()+1;
//        Array<GO> elementList(repeatedMapEpetraPress->NumMyElements());
//        for (unsigned i=0; i<elementList.size(); i++) {
//            elementList[i] = repeatedMapEpetraPress->GID(i) + offsetVelocityMap;
//        }
//        RCP<Map<LO,GO,NO> > repeatedMapPress = MapFactory<LO,GO,NO>::Build(xpetraLib,-1,elementList,0,Comm);

        ArrayRCP<RCP<const Map<LO,GO,NO> > > repeatedMapsVectorConst(2);
        ArrayRCP<RCP<Map<LO,GO,NO> > > repeatedMapsVector(2);

        repeatedMapsVectorConst[0] = repeatedMapVelo.getConst();
        repeatedMapsVectorConst[1] = repeatedMapPress.getConst();

        RCP<Map<LO,GO,NO> > repeatedMap = MergeMapsNonConst(repeatedMapsVectorConst);

        repeatedMapsVector[0] = repeatedMapVelo;
        repeatedMapsVector[1] = repeatedMapPress;


        ////////////////
        // Unique Map //
        ////////////////
        RCP<Map<LO,GO,NO> > uniqueMap = rcp_const_cast<Map<LO,GO,NO> >(BuildUniqueMap<LO,GO,NO>(repeatedMap));


        /////////
        // RHS //
        /////////
        string groupNameRHS = "RHS";
        Epetra_MultiVector *rhsEpetra;
        hDF5IO->Read(groupNameRHS,rhsEpetra);

        RCP<MultiVector<SC,LO,GO,NO> > rhsTmp = ConvertToXpetra<SC,LO,GO,NO>::ConvertMultiVector(xpetraLib,*rhsEpetra,Comm);
        RCP<MultiVector<SC,LO,GO,NO> > rhs = MultiVectorFactory<SC,LO,GO,NO>::Build(uniqueMap,1);
        RCP<Import<LO,GO,NO> > scatter = ImportFactory<LO,GO,NO>::Build(rhsTmp->getMap(),uniqueMap);
        rhs->doImport(*rhsTmp,*scatter,ADD);


        ////////////
        // Matrix //
        ////////////
        string groupNameMatrix = "Matrix";
        Epetra_CrsMatrix *matrixEpetra;
        hDF5IO->Read(groupNameMatrix,matrixEpetra);
        RCP<Matrix<SC,LO,GO,NO> > matrixTmp = ConvertToXpetra<SC,LO,GO,NO>::ConvertMatrix(xpetraLib,*matrixEpetra,Comm);
        RCP<Matrix<SC,LO,GO,NO> > matrix = MatrixFactory<SC,LO,GO,NO>::Build(uniqueMap,matrixTmp->getGlobalMaxNumRowEntries());

        matrix->doImport(*matrixTmp,*scatter,ADD);
        matrix->fillComplete();


        //////////////////////
        // Convert to Thyra //
        //////////////////////
        RCP<MultiVector<SC,LO,GO,NO> > xSolution = MultiVectorFactory<SC,LO,GO,NO>::Build(uniqueMap,1);
        RCP<MultiVector<SC,LO,GO,NO> > xRightHandSide = MultiVectorFactory<SC,LO,GO,NO>::Build(uniqueMap,1);

        xSolution->putScalar(ScalarTraits<SC>::zero());
        xRightHandSide->putScalar(ScalarTraits<SC>::one());

        RCP<CrsMatrixWrap<SC,LO,GO,NO> > tmpCrsWrap = rcp_dynamic_cast<CrsMatrixWrap<SC,LO,GO,NO> >(matrix);

        RCP<const Thyra::LinearOpBase<SC> > K_thyra = ThyraUtils<SC,LO,GO,NO>::toThyra(tmpCrsWrap->getCrsMatrix());

        RCP<Thyra::MultiVectorBase<SC> >thyraX =
        rcp_const_cast<Thyra::MultiVectorBase<SC> >(ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(xSolution));
        RCP<const Thyra::MultiVectorBase<SC> >thyraB = ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(xRightHandSide);


        ///////////////////
        // ParameterList //
        ///////////////////
        RCP<ParameterList> plList = sublist(parameterList,"Preconditioner Types");
        sublist(plList,"FROSch")->set("Overlap",Overlap);

        sublist(plList,"FROSch")->set("Repeated Map Vector",repeatedMapsVector);

        ArrayRCP<DofOrdering> dofOrderings(2);
        ArrayRCP<UN> dofsPerNodeVector(2);
        dofOrderings[0] = NodeWise;
        dofOrderings[1] = NodeWise;
        dofsPerNodeVector[0] = Dimension;
        dofsPerNodeVector[1] = 1;

        sublist(plList,"FROSch")->set("DofOrdering Vector",dofOrderings);
        sublist(plList,"FROSch")->set("DofsPerNode Vector",dofsPerNodeVector);

        Comm->barrier();
        if (Comm->getRank()==0) {
            cout << "##################\n# Parameter List #\n##################" << endl;
            parameterList->print(cout);
            cout << endl;
        }

        Comm->barrier(); if (Comm->getRank()==0) cout << "###################################\n# Stratimikos LinearSolverBuilder #\n###################################\n" << endl;
        Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
        Stratimikos::enableFROSch<LO,GO,NO>(linearSolverBuilder);
        linearSolverBuilder.setParameterList(parameterList);

        Comm->barrier(); if (Comm->getRank()==0) cout << "######################\n# Thyra PrepForSolve #\n######################\n" << endl;

        RCP<LinearOpWithSolveFactoryBase<SC> > lowsFactory =
        linearSolverBuilder.createLinearSolveStrategy("");

        lowsFactory->setOStream(out);
        lowsFactory->setVerbLevel(VERB_HIGH);

        Comm->barrier(); if (Comm->getRank()==0) cout << "###########################\n# Thyra LinearOpWithSolve #\n###########################" << endl;

        RCP<LinearOpWithSolveBase<SC> > lows =
        linearOpWithSolve(*lowsFactory, K_thyra);

        Comm->barrier(); if (Comm->getRank()==0) cout << "\n#########\n# Solve #\n#########" << endl;
        SolveStatus<double> status =
        solve<double>(*lows, Thyra::NOTRANS, *thyraB, thyraX.ptr());

        Comm->barrier(); if (Comm->getRank()==0) cout << "\n#############\n# Finished! #\n#############" << endl;
    }
#endif

    CommWorld->barrier();
    stackedTimer->stop("Thyra Stokes Test");
    StackedTimer::OutputOptions options;
    options.output_fraction = options.output_histogram = options.output_minmax = true;
    stackedTimer->report(*out,CommWorld,options);

    return(EXIT_SUCCESS);

}



