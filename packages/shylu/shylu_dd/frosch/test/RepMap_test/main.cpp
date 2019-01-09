#include <mpi.h>
#include <Epetra_MpiComm.h>

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Import.h"
#include "Epetra_SerialDenseMatrix.h"

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
#include "FROSch_Tools_def.hpp"

typedef unsigned UN;
typedef unsigned UN;
typedef double SC;
typedef int LO;
typedef int GO;
typedef Kokkos::Compat::KokkosSerialWrapperNode NO;
//typedef KokkosClassic::DefaultNode::DefaultNodeType  EpetraNode;
//typedef EpetraNode NO;

using namespace std;
using namespace Teuchos;
using namespace Xpetra;
using namespace FROSch;
using namespace Belos;

int main(int argc, char *argv[])
{
    MPI_Init(&argc,&argv);
    
    {
        int NumMyElements = 0;                                          // NODES assigned to this processor
        int NumMyExternalElements = 0;                             // nodes used by this proc, but not hosted
        int NumMyTotalElements = 0;
        int FE_NumMyElements = 0;                                    // TRIANGLES assigned to this processor
        Teuchos::Array<GO> MyGlobalElements(8);           // nodes assigned to this processor
        Teuchos::Array<GO> FEelements(8);                     // elements assigned to this proc
        Epetra_IntSerialDenseMatrix E;                             // store the element graph connectivity
        Epetra_IntSerialDenseMatrix ElementNodeList;
        Teuchos::Array<Teuchos::Array<GO> > NodesInElement(8);
        
        Teuchos::RCP<const Teuchos::Comm<int> > TeuchosComm = Teuchos::rcp(new Teuchos::MpiComm<int> (MPI_COMM_WORLD));
        
        
        Teuchos::RCP<Teuchos::FancyOStream> fancy = fancyOStream(Teuchos::rcpFromRef(std::cout));
        
        
        int MyPID=TeuchosComm->getRank();
        
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
        TeuchosComm->barrier();TeuchosComm->barrier();TeuchosComm->barrier();
        if(MyPID == 0) cout<<"Set Up Map finished\n";
        
        Teuchos::RCP<Xpetra::Map<LO,GO,NO> > Mapg = Xpetra::MapFactory<LO,GO,NO>::Build(Xpetra::UseTpetra,-1,FEelements(),0,TeuchosComm);
        
        Teuchos::RCP<Xpetra::TpetraMap<LO,GO,NO> > MapT = Teuchos::rcp(new Xpetra::TpetraMap<LO,GO,NO>  (32,FEelements(),0,TeuchosComm));
        
        TeuchosComm->barrier();TeuchosComm->barrier();TeuchosComm->barrier();
        if(MyPID == 0) cout<<" Map build\n";
        
        
        switch( MyPID ) {
                case 0:
                E.Shape(FE_NumMyElements,4);
                //graph connectivity
                E(0,0) = 1; E(0,1) = 9;   E(0,2) = 0;  E(0,3) = 0;
                E(1,0) = 0; E(1,1) = 2;   E(1,2) = 1;  E(1,3) = 1;
                E(2,0) = 1; E(2,1) = 3;   E(2,2) = 11; E(2,3) = 2;
                E(3,0) = 2; E(3,1) = 4;   E(3,2) = 3;  E(3,3) = 3;
                E(4,0) = 3; E(4,1) = 5;   E(4,2) = 13; E(4,3) = 4;
                E(5,0) = 4; E(5,1) = 6;   E(5,2) = 5;  E(5,3) = 5;
                E(6,0) = 5; E(6,1) = 7;   E(6,2) = 15; E(6,3) = 6;
                E(7,0) = 6; E(7,1) = 7;   E(7,2) = 7;  E(7,3) = 7;
                break;
                case 1:
                E.Shape(FE_NumMyElements,4);
                
                
                //graph connectivity
                E(0,0) = 9;   E(0,1) = 17;  E(0,2) = 8;    E(0,3) = 8;
                E(1,0) = 0;   E(1,1) = 8;   E(1,2) = 10;   E(1,3) = 9;
                E(2,0) = 9;   E(2,1) = 11;  E(2,2) = 19;   E(2,3) = 10;
                E(3,0) = 2;   E(3,1) = 10;  E(3,2) = 12;   E(3,3) = 11;
                E(4,0) = 11;  E(4,1) = 13;  E(4,2) = 21;   E(4,3) = 12;
                E(5,0) = 4;   E(5,1) = 12;  E(5,2) = 14;   E(5,3) = 13;
                E(6,0) = 13;  E(6,1) = 15;  E(6,2) = 23;   E(6,3) = 14;
                E(7,0) = 6;   E(7,1) = 14;  E(7,2) = 15;   E(7,3) = 15;
                break;
                case 2:
                E.Shape(FE_NumMyElements,4);
                
                //graph connectivity
                E(0,0) = 17; E(0,1) = 25; E(0,2) = 16;   E(0,3) = 16;
                E(1,0) = 8;  E(1,1) = 16; E(1,2) = 18;   E(1,3) = 17;
                E(2,0) = 17; E(2,1) = 19; E(2,2) = 27;   E(2,3) = 18;
                E(3,0) = 10; E(3,1) = 18; E(3,2) = 20;   E(3,3) = 19;
                E(4,0) = 19; E(4,1) = 21; E(4,2) = 29;   E(4,3) = 20;
                E(5,0) = 12; E(5,1) = 20; E(5,2) = 22;   E(5,3) = 21;
                E(6,0) = 21; E(6,1) = 23; E(6,2) = 31;   E(6,3) = 22;
                E(7,0) = 14; E(7,1) = 22; E(7,2) = 23;   E(7,3) = 23;
                break;
                case 3:
                E.Shape(FE_NumMyElements,4);
                //graph connectivity
                E(0,0) = 25; E(0,1) = 24; E(0,2) = 24;   E(0,3) = 24;
                E(1,0) = 16; E(1,1) = 24;  E(1,2) = 26;  E(1,3) = 25;
                E(2,0) = 25; E(2,1) = 27;  E(2,2) = 26;  E(2,3) = 26;
                E(3,0) = 18; E(3,1) = 26;  E(3,2) = 28;  E(3,3) = 27;
                E(4,0) = 27; E(4,1) = 29;  E(4,2) = 28;  E(4,3) = 28;
                E(5,0) = 20; E(5,1) = 28;  E(5,2) = 30;  E(5,3) = 29;
                E(6,0) = 29; E(6,1) = 31;  E(6,2) = 30;  E(6,3) = 30;
                E(7,0) = 22; E(7,1) = 30;  E(7,2) = 31;  E(7,3) = 31;
                break;
        }
        TeuchosComm->barrier();TeuchosComm->barrier();TeuchosComm->barrier();
        if(MyPID == 0) cout<<"Element Connectivity\n";
        
        
        const Teuchos::RCP<Xpetra::CrsGraph<LO,GO,NO> > Xgraph = Xpetra::CrsGraphFactory<LO,GO,NO>::Build(Mapg,4);
        const int MaxNnzRow = 8;
        const int MaxNnzRowg = 4 ;
        
        int Element, MyRow, GlobalRow, GlobalCol;
        Teuchos::ArrayView<int> array;
        Teuchos::Array<int> v(1);
        v[0] = 1.0;
        Teuchos::Array<int> gcol(1);
        Teuchos::Array<int> indi(3);
        int numGentry = 4;
        
        TeuchosComm->barrier();TeuchosComm->barrier();TeuchosComm->barrier();
        if(MyPID == 0) cout<<"329\n";
        
        for(  int Element=0 ; Element<FE_NumMyElements ; Element++ )
        {
            // variables used inside
            int GlobalRow;
            int MyRow;
            int GlobalCol;
            int count;
            Teuchos::Array<int> indices(4);
            // get the global and local number of this row
            MyRow = Element;
            GlobalRow = MyRow+MyPID*8;
            
            if( MyRow < FE_NumMyElements) {
                for( int j=0 ; j<4 ; j++ ) {
                    gcol[0] = E(Element,j);
                    indices[j] = E(Element,j);
                }
                numGentry = indices.size();
                Xgraph->insertGlobalIndices(GlobalRow,indices());
            }
        }
        Xgraph->fillComplete();
        TeuchosComm->barrier();TeuchosComm->barrier();TeuchosComm->barrier();
        if(MyPID == 0) cout<<"Xgraph\n";
        
        
        
        
        Teuchos::RCP<const Xpetra::Map<LO, GO, NO> > map = Xpetra::MapFactory<LO, GO, NO>::createUniformContigMap(Xpetra::UseTpetra,Mapg->getGlobalNumElements(), TeuchosComm);
        
        
        //Node ElementList in B---------------------------------------------------------------------
        Teuchos::RCP<Xpetra::TpetraCrsMatrix<GO> > B = Teuchos::rcp(new Xpetra::TpetraCrsMatrix<GO> (Mapg, 3));
        
        const size_t numMyElements = MapT->getNodeNumElements ();
        size_t numLocalElements = MapT->getMaxLocalIndex();
        Teuchos::ArrayView<const GO> myGlobalElements = MapT->getNodeElementList();
        
        std::vector<GO> col_vec(3);
        std::vector<int> val_ele(3);
        for(int i = 0;i<3;i++)
        {
            col_vec.at(i) =i;
            
        }
        Teuchos::ArrayView<GO> cols(col_vec);
        
        for (size_t i = 0; i < numMyElements; ++i){
            for(int j = 0 ;j<3;j++){
                val_ele.at(j) = NodesInElement.at(i).at(j);
            }
            Teuchos::ArrayView<GO> vals(val_ele);
            //NodeEleList->insertGlobalValues(i*(TeuchosComm->getRank()+1),cols,vals);
            B->insertGlobalValues (myGlobalElements[i],cols, vals);
        }
        B->fillComplete();
        //B->describe(*fancy,Teuchos::VERB_EXTREME);
        string xmlFile = "xpetra_ParameterList.xml";
        Teuchos::RCP<Teuchos::ParameterList> parameterList = Teuchos::getParametersFromXmlFile(xmlFile);
        
        Teuchos::RCP<Xpetra::Map<LO,GO,NO> > RepeatedMap = FROSch::BuildRepMap_Zoltan<SC,LO,GO,NO>(Xgraph,B,parameterList,TeuchosComm);
        //RepeatedMap->describe(*fancy,Teuchos::VERB_EXTREME);
        const Teuchos::RCP<const Xpetra::Map<LO,GO,NO> > cRep = RepeatedMap;
        Teuchos::RCP<Xpetra::Map<LO,GO,NO> > UniqueMap = FROSch::BuildUniqueMap(cRep);
        
    }
    
    MPI_Finalize();
    
    return(EXIT_SUCCESS);

}

void compute_loc_matrix( double *x_triangle, double *y_triangle,
                        Epetra_SerialDenseMatrix & Ke )
{
    int    ii, jj;
    double det_J;
    double xa, ya, xb, yb, xc, yc;
    xa = x_triangle[0];
    xb = x_triangle[1];
    xc = x_triangle[2];
    ya = y_triangle[0];
    yb = y_triangle[1];
    yc = y_triangle[2];
    Ke(0,0) = (yc-yb)*(yc-yb) + (xc-xb)*(xc-xb);
    Ke(0,1) = (yc-yb)*(ya-yc) + (xc-xb)*(xa-xc);
    Ke(0,2) = (yb-ya)*(yc-yb) + (xb-xa)*(xc-xb);
    Ke(1,0) = (yc-yb)*(ya-yc) + (xc-xb)*(xa-xc);
    Ke(1,1) = (yc-ya)*(yc-ya) + (xc-xa)*(xc-xa);
    Ke(1,2) = (ya-yc)*(yb-ya) + (xa-xc)*(xb-xa);
    Ke(2,0) = (yb-ya)*(yc-yb) + (xb-xa)*(xc-xb);
    Ke(2,1) = (ya-yc)*(yb-ya) + (xa-xc)*(xb-xa);
    Ke(2,2) = (yb-ya)*(yb-ya) + (xb-xa)*(xb-xa);
    det_J = (xb-xa)*(yc-ya)-(xc-xa)*(yb-ya);
    det_J = 2*det_J;
    if( det_J<0 ) det_J *=-1;
    for (ii = 0; ii < 3; ii++) {
        for (jj = 0; jj < 3; jj++) {
            Ke(ii,jj) = Ke(ii,jj) / det_J;
        }
    }
    return;
} /* compute_loc_matrix */




