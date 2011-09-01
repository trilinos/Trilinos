#include"Isorropia_EpetraMatcher.hpp"
#include"Isorropia_EpetraRedistributor.hpp"

#ifdef HAVE_EPETRAEXT
#include "EpetraExt_Reindex_CrsMatrix.h" 
#include "EpetraExt_CrsMatrixIn.h"
#endif


using namespace std;

int main(int argc, char** argv) {

#ifdef ISORROPIA_HAVE_OMP
	if(argc>2)
	{	
		int rc=0;
		int localProc = 0;
		
		#ifdef HAVE_EPETRAEXT
		/*#ifdef HAVE_MPI
			int numProcs;
		 	MPI_Init(&argc, &argv);
		  	MPI_Comm_rank(MPI_COMM_WORLD, &localProc);
		  	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
		  	const Epetra_MpiComm Comm(MPI_COMM_WORLD);
		  	const Epetra_MpiComm Comm;
		#else*/
		  const Epetra_SerialComm Comm;
		//#endif

		  Epetra_CrsMatrix *matrixPtr;
		  rc = EpetraExt::MatrixMarketFileToCrsMatrix(argv[1], Comm, matrixPtr);
		  
		  if (rc < 0){
			 if (localProc==0){
				cout << "error reading input file" << std::endl << "FAIL" << std::endl;
			 }
			 exit(1);
		  }
		#else
		  fail = 0;
		  if (localProc == 0){
			 cout << "Test not run because it requires EPETRA_EXT" << std::endl;
		  }
		#endif
		
		Teuchos::ParameterList paramlist;
		paramlist.set("Matching Algorithm",argv[2]);
        Isorropia::Epetra::Isorropia_EpetraMatcher pm(matrixPtr,paramlist);
		
        //Teuchos::RCP<const Epetra_CrsMatrix> r(Teuchos::RCP<const
        //Epetra_CrsMatrix>(matrixPtr,true));
        //Isorropia::Epetra::Isorropia_EpetraMatcher pm(r,paramlist);
        
        pm.match();
        
        /*std::cout<<*matrixPtr<<std::endl;
        Epetra_Map * map = pm.getPermutedRowMap();
        Teuchos::RCP<Epetra_Map> rcpMap(Teuchos::RCP<
        Epetra_Map>(map,true));
        Isorropia::Epetra::Redistributor redist(rcpMap);
        Teuchos::RCP<Epetra_CrsMatrix> myMat=redist.redistribute(*matrixPtr);
        std::cout<<*myMat<<std::endl;*/

        //Epetra_Map defMap2(-1, 5, 0,matrixPtr->Comm());
        //EpetraExt::ViewTransform<Epetra_CrsMatrix> * ReIdx_MatTrans2 =
                               //new EpetraExt::CrsMatrix_Reindex( defMap2 ); 
        //Epetra_CrsMatrix t2S = (*ReIdx_MatTrans2)(*myMat);
        //ReIdx_MatTrans2->fwd(); 
        //std::cout<<t2S<<std::endl;
	}
	else
    {
		cout<<endl<<" Usage: ./Isorropia_parallel_matching.exe <mtx file>" <<
             " <Algorithm>" << endl;
        cout << "\t Algorithm: PHK, PHKDW, PDFS,PPF" << endl << endl;
    }
#else
    cout << "Matching in Isorropia requires OpenMP." << endl;
    cout << "Please recompile with OpenMP enabled and try again" << endl;
#endif
	
	return 0;
}
