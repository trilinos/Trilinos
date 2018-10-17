#include <mpi.h>
#include <Epetra_MpiComm.h>

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>

#include <Galeri_Maps.h>
#include <Galeri_CrsMatrices.h>
#include "Galeri_Utils.h"


#include <Teuchos_CommandLineProcessor.hpp>
#include <Teuchos_XMLParameterListCoreHelpers.hpp>

#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_EpetraCrsMatrix.hpp>

#include <BelosOperatorT.hpp>
#include <BelosXpetraAdapter.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>

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
#include "Stratimikos_FROSchXpetra.hpp"

// Xpetra include
#include <Xpetra_Parameters.hpp>

// FROSCH thyra includes
#include "Thyra_FROSchLinearOp_def.hpp"
#include "Thyra_FROSchXpetraFactory_def.hpp"
#include "FROSch_Tools_def.hpp"

typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef int GO;
typedef Kokkos::Compat::KokkosSerialWrapperNode EpetraNode; // Hier Default verwenden???
typedef EpetraNode NO;

using namespace std;
using namespace Teuchos;
using namespace Xpetra;
using namespace FROSch;
using namespace Belos;

int main(int argc, char *argv[])
{
    MPI_Init(&argc,&argv);
    
    {
          Epetra_MpiComm Comm(MPI_COMM_WORLD);
        int NumMyElements = 0;         // NODES assigned to this processor
        int NumMyExternalElements = 0; // nodes used by this proc, but not hosted
        int NumMyTotalElements = 0;
        int FE_NumMyElements = 0;      // TRIANGLES assigned to this processor
        Teuchos::Array<GO> MyGlobalElements(8);    // nodes assigned to this processor
        Teuchos::Array<GO> FEelements(8);
        // elements assigned to this proc
        Teuchos::Array<Teuchos::Array<long long> > E(8); // store the element graph connectivity
        Epetra_IntSerialDenseMatrix ElementNodeList;
        Teuchos::Array<Teuchos::Array<long long> > NodesInElement(8);
        
        
        Teuchos::RCP<const Teuchos::Comm<int> > TeuchosComm = Teuchos::rcp(new Teuchos::MpiComm<int> (MPI_COMM_WORLD));
        
        
        string xmlFile = "xpetra_ParameterList.xml";
        Teuchos::RCP<Teuchos::ParameterList> parameterList = Teuchos::getParametersFromXmlFile(xmlFile);
        
        Teuchos::RCP<Teuchos::FancyOStream> fancy = fancyOStream(Teuchos::rcpFromRef(std::cout));
        
        
        int MyPID=Comm.MyPID();
        switch (MyPID){
            case 0:
                NumMyElements = 6;
                NumMyExternalElements = 4;
                NumMyTotalElements = NumMyElements + NumMyExternalElements;
                FE_NumMyElements = 8;
                
                
                
                MyGlobalElements[0] = 0;
                MyGlobalElements[1] = 1;
                MyGlobalElements[2] = 2;
                MyGlobalElements[3] = 3;
                MyGlobalElements[4] = 4;
                MyGlobalElements[5] = 5;
                MyGlobalElements[6] = 6;
                MyGlobalElements[7] = 7;
                MyGlobalElements[8] = 8;
                MyGlobalElements[9] = 9;
                
                FEelements[0] = 0;
                FEelements[1] = 1;
                FEelements[2] = 2;
                FEelements[3] = 3;
                FEelements[4] = 4;
                FEelements[5] = 5;
                FEelements[6] = 6;
                FEelements[7] = 7;
                
                NodesInElement.at(0).resize(3);
                NodesInElement.at(1).resize(3);
                NodesInElement.at(2).resize(3);
                NodesInElement.at(3).resize(3);
                NodesInElement.at(4).resize(3);
                NodesInElement.at(5).resize(3);
                NodesInElement.at(6).resize(3);
                NodesInElement.at(7).resize(3);
                
                
                NodesInElement.at(0).at(0) = 0 ; NodesInElement.at(0).at(1) = 5; NodesInElement.at(0).at(2) = 6;
                
                NodesInElement.at(1).at(0) = 0 ; NodesInElement.at(1).at(1) = 1; NodesInElement.at(1).at(2) = 6;
                NodesInElement.at(2).at(0) = 1 ; NodesInElement.at(2).at(1) = 6; NodesInElement.at(2).at(2) = 7;
                NodesInElement.at(3).at(0) = 1 ; NodesInElement.at(3).at(1) = 2; NodesInElement.at(3).at(2) = 7;
                NodesInElement.at(4).at(0) = 2 ; NodesInElement.at(4).at(1) = 7; NodesInElement.at(4).at(2) = 8;
                NodesInElement.at(5).at(0) = 2 ; NodesInElement.at(5).at(1) = 3; NodesInElement.at(5).at(2) = 8;
                NodesInElement.at(6).at(0) = 3 ; NodesInElement.at(6).at(1) = 8; NodesInElement.at(6).at(2) = 9;
                NodesInElement.at(7).at(0) = 3 ; NodesInElement.at(7).at(1) = 4; NodesInElement.at(7).at(2) = 9;
                break;
                
            case 1:
                NumMyElements = 6;
                NumMyExternalElements = 4;
                NumMyTotalElements = NumMyElements + NumMyExternalElements;
                FE_NumMyElements = 8;
                
                NodesInElement.at(0).resize(3);
                NodesInElement.at(1).resize(3);
                NodesInElement.at(2).resize(3);
                NodesInElement.at(3).resize(3);
                NodesInElement.at(4).resize(3);
                NodesInElement.at(5).resize(3);
                NodesInElement.at(6).resize(3);
                NodesInElement.at(7).resize(3);
                
                
                MyGlobalElements[0] = 6;
                MyGlobalElements[1] = 7;
                MyGlobalElements[2] = 8;
                MyGlobalElements[3] = 9;
                MyGlobalElements[4] = 10;
                MyGlobalElements[5] = 11;
                MyGlobalElements[6] = 5;
                MyGlobalElements[7] = 12;
                MyGlobalElements[8] = 13;
                MyGlobalElements[9] = 14;
                
                FEelements[0] = 8;
                FEelements[1] = 9;
                FEelements[2] = 10;
                FEelements[3] = 11;
                FEelements[4] = 12;
                FEelements[5] = 13;
                FEelements[6] = 14;
                FEelements[7] = 15;
                NodesInElement.at(0).at(0) = 5 ; NodesInElement.at(0).at(1) = 10; NodesInElement.at(0).at(2) = 11;
                NodesInElement.at(1).at(0) = 5 ; NodesInElement.at(1).at(1) = 6; NodesInElement.at(1).at(2) = 11;
                NodesInElement.at(2).at(0) = 6 ; NodesInElement.at(2).at(1) = 11; NodesInElement.at(2).at(2) = 12;
                NodesInElement.at(3).at(0) = 6 ; NodesInElement.at(3).at(1) = 7; NodesInElement.at(3).at(2) = 12;
                NodesInElement.at(4).at(0) = 7 ; NodesInElement.at(4).at(1) = 12; NodesInElement.at(4).at(2) = 13;
                NodesInElement.at(5).at(0) = 7 ; NodesInElement.at(5).at(1) = 8; NodesInElement.at(5).at(2) = 13;
                NodesInElement.at(6).at(0) = 8 ; NodesInElement.at(6).at(1) = 13; NodesInElement.at(6).at(2) = 14;
                NodesInElement.at(7).at(0) = 8 ; NodesInElement.at(7).at(1) = 9; NodesInElement.at(7).at(
                                                                                                         
                                                                                                         2) = 14;
                break;
            case 2:
                NumMyElements = 6;
                NumMyExternalElements = 4;
                NumMyTotalElements = NumMyElements + NumMyExternalElements;
                FE_NumMyElements = 8;
                
                NodesInElement.at(0).resize(3);
                NodesInElement.at(1).resize(3);
                NodesInElement.at(2).resize(3);
                NodesInElement.at(3).resize(3);
                NodesInElement.at(4).resize(3);
                NodesInElement.at(5).resize(3);
                NodesInElement.at(6).resize(3);
                NodesInElement.at(7).resize(3);
                
                MyGlobalElements[0] = 12;
                MyGlobalElements[1] = 13;
                MyGlobalElements[2] = 14;
                MyGlobalElements[3] = 15;
                MyGlobalElements[4] = 16;
                MyGlobalElements[5] = 17;
                MyGlobalElements[6] = 10;
                MyGlobalElements[7] = 11;
                MyGlobalElements[8] = 18;
                MyGlobalElements[9] = 19;
                
                FEelements[0] = 16;
                FEelements[1] = 17;
                FEelements[2] = 18;
                FEelements[3] = 19;
                FEelements[4] = 20;
                FEelements[5] = 21;
                FEelements[6] = 22;
                FEelements[7] = 23;
                
                NodesInElement.at(0).at(0) = 10 ; NodesInElement.at(0).at(1) = 15; NodesInElement.at(0).at(2) = 16;
                NodesInElement.at(1).at(0) = 10 ; NodesInElement.at(1).at(1) = 11; NodesInElement.at(1).at(2) = 16;
                NodesInElement.at(2).at(0) = 11 ; NodesInElement.at(2).at(1) = 16; NodesInElement.at(2).at(2) = 17;
                NodesInElement.at(3).at(0) = 11 ; NodesInElement.at(3).at(1) = 12; NodesInElement.at(3).at(2) = 17;
                NodesInElement.at(4).at(0) = 12 ; NodesInElement.at(4).at(1) = 17; NodesInElement.at(4).at(2) = 18;
                NodesInElement.at(5).at(0) = 12 ; NodesInElement.at(5).at(1) = 13; NodesInElement.at(5).at(2) = 18;
                NodesInElement.at(6).at(0) = 13 ; NodesInElement.at(6).at(1) = 18; NodesInElement.at(6).at(2) = 19;
                NodesInElement.at(7).at(0) = 13 ; NodesInElement.at(7).at(1) = 14; NodesInElement.at(7).at(2) = 19;
                break;
            case 3:
                NumMyElements = 7;
                NumMyExternalElements = 3;
                NumMyTotalElements = NumMyElements + NumMyExternalElements;
                FE_NumMyElements = 8;
                
                NodesInElement.at(0).resize(3);
                NodesInElement.at(1).resize(3);
                NodesInElement.at(2).resize(3);
                NodesInElement.at(3).resize(3);
                NodesInElement.at(4).resize(3);
                NodesInElement.at(5).resize(3);
                NodesInElement.at(6).resize(3);
                NodesInElement.at(7).resize(3);
                
                
                
                MyGlobalElements[0] = 18;
                MyGlobalElements[1] = 19;
                MyGlobalElements[2] = 20;
                MyGlobalElements[3] = 21;
                MyGlobalElements[4] = 22;
                MyGlobalElements[5] = 23;
                MyGlobalElements[6] = 24;
                MyGlobalElements[7] = 15;
                MyGlobalElements[8] = 16;
                MyGlobalElements[9] = 17;
                
                FEelements[0] = 24;
                FEelements[1] = 25;
                FEelements[2] = 26;
                FEelements[3] = 27;
                FEelements[4] = 28;
                FEelements[5] = 29;
                FEelements[6] = 30;
                FEelements[7] = 31;
                
                NodesInElement.at(0).at(0) = 15 ; NodesInElement.at(0).at(1) = 20; NodesInElement.at(0).at(2) = 21;
                NodesInElement.at(1).at(0) = 15 ; NodesInElement.at(1).at(1) = 16; NodesInElement.at(1).at(2) = 21;
                NodesInElement.at(2).at(0) = 16 ; NodesInElement.at(2).at(1) = 21; NodesInElement.at(2).at(2) = 22;
                NodesInElement.at(3).at(0) = 16 ; NodesInElement.at(3).at(1) = 17; NodesInElement.at(3).at(2) = 22;
                NodesInElement.at(4).at(0) = 17 ; NodesInElement.at(4).at(1) = 22; NodesInElement.at(4).at(2) = 23;
                NodesInElement.at(5).at(0) = 17 ; NodesInElement.at(5).at(1) = 18; NodesInElement.at(5).at(2) = 23;
                NodesInElement.at(6).at(0) = 18 ; NodesInElement.at(6).at(1) = 23; NodesInElement.at(6).at(2) = 24;
                NodesInElement.at(7).at(0) = 18 ; NodesInElement.at(7).at(1) = 19; NodesInElement.at(7).at(2) = 24;
                
                break;
        }
        Comm.Barrier();    Comm.Barrier();    Comm.Barrier();
        if(MyPID == 0) cout<<"Set Up Map finished\n";
        
        Teuchos::RCP<Xpetra::Map<LO,GO,NO> > Mapg = Xpetra::MapFactory<LO,GO,NO>::Build(Xpetra::UseEpetra,-1,FEelements(),0,TeuchosComm);
        
        
        Teuchos::RCP<Xpetra::MultiVector<GO,LO,GO,NO> >
        NodeEleList  = Xpetra::MultiVectorFactory<GO,LO,GO,NO>::Build(Mapg,3);
        for(int i = 0;i<8;i++){
            for(int j = 0;j<3;j++)
                NodeEleList->replaceLocalValue(i,j,NodesInElement.at(i).at(j));
        }
    
    switch( MyPID ) {
        case 0:
            E.at(0).resize(4);
            E.at(1).resize(4);
            E.at(2).resize(4);
            E.at(3).resize(4);
            E.at(4).resize(4);
            E.at(5).resize(4);
            E.at(6).resize(4);
            E.at(7).resize(4);
            //graph connectivity
            E.at(0).at(0) = 1; E.at(0).at(1) = 9;   E.at(0).at(2) = 0;  E.at(0).at(3) = 0;
            E.at(1).at(0) = 0; E.at(1).at(1) = 2;   E.at(1).at(2) = 1;  E.at(1).at(3) = 1;
            E.at(2).at(0) = 1; E.at(2).at(1) = 3;   E.at(2).at(2) = 11; E.at(2).at(3) = 2;
            E.at(3).at(0) = 2; E.at(3).at(1) = 4;   E.at(3).at(2) = 3;  E.at(3).at(3) = 3;
            E.at(4).at(0) = 3; E.at(4).at(1) = 5;   E.at(4).at(2) = 13; E.at(4).at(3) = 4;
            E.at(5).at(0) = 4; E.at(5).at(1) = 6;   E.at(5).at(2) = 5;  E.at(5).at(3) = 5;
            E.at(6).at(0) = 5; E.at(6).at(1) = 7;   E.at(6).at(2) = 15; E.at(6).at(3) = 6;
            E.at(7).at(0) = 6; E.at(7).at(1) = 7;   E.at(7).at(2) = 7;  E.at(7).at(3) = 7;
            break;
        case 1:
            E.at(0).resize(4);
            E.at(1).resize(4);
            E.at(2).resize(4);
            E.at(3).resize(4);
            E.at(4).resize(4);
            E.at(5).resize(4);
            E.at(6).resize(4);
            E.at(7).resize(4);
            
            //graph connectivity
            E.at(0).at(0) = 9;   E.at(0).at(1) = 17;  E.at(0).at(2) = 8;    E.at(0).at(3) = 8;
            E.at(1).at(0) = 0;   E.at(1).at(1) = 8;   E.at(1).at(2) = 10;   E.at(1).at(3) = 9;
            E.at(2).at(0) = 9;   E.at(2).at(1) = 11;  E.at(2).at(2) = 19;   E.at(2).at(3) = 10;
            E.at(3).at(0) = 2;   E.at(3).at(1) = 10;  E.at(3).at(2) = 12;   E.at(3).at(3) = 11;
            E.at(4).at(0) = 11;  E.at(4).at(1) = 13;  E.at(4).at(2) = 21;   E.at(4).at(3) = 12;
            E.at(5).at(0) = 4;   E.at(5).at(1) = 12;  E.at(5).at(2) = 14;   E.at(5).at(3) = 13;
            E.at(6).at(0) = 13;  E.at(6).at(1) = 15;  E.at(6).at(2) = 23;   E.at(6).at(3) = 14;
            E.at(7).at(0) = 6;   E.at(7).at(1) = 14;  E.at(7).at(2) = 15;   E.at(7).at(3) = 15;
            break;
        case 2:
            E.at(0).resize(4);
            E.at(1).resize(4);
            E.at(2).resize(4);
            E.at(3).resize(4);
            E.at(4).resize(4);
            E.at(5).resize(4);
            E.at(6).resize(4);
            E.at(7).resize(4);
            
            //graph connectivity
            E.at(0).at(0) = 17; E.at(0).at(1) = 25; E.at(0).at(2) = 16;   E.at(0).at(3) = 16;
            E.at(1).at(0) = 8;  E.at(1).at(1) = 16; E.at(1).at(2) = 18;   E.at(1).at(3) = 17;
            E.at(2).at(0) = 17; E.at(2).at(1) = 19; E.at(2).at(2) = 27;   E.at(2).at(3) = 18;
            E.at(3).at(0) = 10; E.at(3).at(1) = 18; E.at(3).at(2) = 20;   E.at(3).at(3) = 19;
            E.at(4).at(0) = 19; E.at(4).at(1) = 21; E.at(4).at(2) = 29;   E.at(4).at(3) = 20;
            E.at(5).at(0) = 12; E.at(5).at(1) = 20; E.at(5).at(2) = 22;   E.at(5).at(3) = 21;
            E.at(6).at(0) = 21; E.at(6).at(1) = 23; E.at(6).at(2) = 31;   E.at(6).at(3) = 22;
            E.at(7).at(0) = 14; E.at(7).at(1) = 22; E.at(7).at(2) = 23;   E.at(7).at(3) = 23;
            break;
        case 3:
            E.at(0).resize(4);
            E.at(1).resize(4);
            E.at(2).resize(4);
            E.at(3).resize(4);
            E.at(4).resize(4);
            E.at(5).resize(4);
            E.at(6).resize(4);
            E.at(7).resize(4);

            //graph connectivity
            E.at(0).at(0) = 25; E.at(0).at(1) = 24;  E.at(0).at(2) = 24;  E.at(0).at(3) = 24;
            E.at(1).at(0) = 16; E.at(1).at(1) = 24;  E.at(1).at(2) = 26;  E.at(1).at(3) = 25;
            E.at(2).at(0) = 25; E.at(2).at(1) = 27;  E.at(2).at(2) = 26;  E.at(2).at(3) = 26;
            E.at(3).at(0) = 18; E.at(3).at(1) = 26;  E.at(3).at(2) = 28;  E.at(3).at(3) = 27;
            E.at(4).at(0) = 27; E.at(4).at(1) = 29;  E.at(4).at(2) = 28;  E.at(4).at(3) = 28;
            E.at(5).at(0) = 20; E.at(5).at(1) = 28;  E.at(5).at(2) = 30;  E.at(5).at(3) = 29;
            E.at(6).at(0) = 29; E.at(6).at(1) = 31;  E.at(6).at(2) = 30;  E.at(6).at(3) = 30;
            E.at(7).at(0) = 22; E.at(7).at(1) = 30;  E.at(7).at(2) = 31;  E.at(7).at(3) = 31;
            break;
    }
      
        Teuchos::RCP<Xpetra::Matrix<GO,LO,GO,NO> >
        connection  = Xpetra::MatrixFactory<GO,LO,GO,NO>::Build(Mapg,4);
        for(int i = 0;i<8;i++){
            for(int j = 0;j<3;j++)
                connection->replaceLocalValues(i,j,E.at(i).at(j));
        }
        
        Teuchos::RCP<Xpetra::Map<LO,GO,NO> > RepeatedMap = FROSch::BuildRepMap_Zoltan<SC,LO,GO,NO>
        (connection,NodeEleList,parameterList);
        
    
    
        
    }
    
    MPI_Finalize();
    
    return(EXIT_SUCCESS);

}



