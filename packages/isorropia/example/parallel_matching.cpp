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
#ifdef HAVE_EPETRAEXT
		  const Epetra_SerialComm Comm;

		  Epetra_CrsMatrix *matrixPtr;
		  rc = EpetraExt::MatrixMarketFileToCrsMatrix(argv[1], Comm, matrixPtr);
		  
		  if (rc < 0){
            cout << "error reading input file" << std::endl;
            return 1;
		  }
		
		Teuchos::ParameterList paramlist;
		paramlist.set("Matching Algorithm",argv[2]);
        Isorropia::Epetra::Isorropia_EpetraMatcher pm(matrixPtr,paramlist);
		
        //Teuchos::RCP<const Epetra_CrsMatrix> r(Teuchos::RCP<const
        //Epetra_CrsMatrix>(matrixPtr,true));
        //Isorropia::Epetra::Isorropia_EpetraMatcher pm(r,paramlist);
        
        pm.match();

        cout << endl << "Original Matrix:" << endl;
        std::cout<<*matrixPtr<<std::endl;

        Teuchos::RCP<Epetra_CrsMatrix> perm_matrix =
                                        pm.applyColumnPermutation();
        cout << endl << "After Column permutation:" << endl;
        cout << *perm_matrix << endl;

        perm_matrix = pm.applyRowPermutation();
        cout << endl << "After Row permutation:" << endl;
        cout << *perm_matrix << endl;

#else
		 fail = 0;
         cout << "Matching test requires EpetraExt" << std::endl;
         return 1;
#endif
	}
	else
    {
		cout<<endl<<" Usage: ./Isorropia_parallel_matching.exe <mtx file>" <<
             " <Algorithm>" << endl;
        cout << "\t Algorithm: PHK, PHKDW, PDFS,PPF" << endl;
        cout << "Requires Isorropia to be compiled with OpenMP, OMP_NUM_THREADS"
              << "set to at least one and the test requires EpetraExt." <<
               endl << endl;
    }
#else
    cout << "Matching in Isorropia requires OpenMP." << endl;
    cout << "Please recompile with OpenMP enabled and try again" << endl;
#endif
	
	return 0;
}
