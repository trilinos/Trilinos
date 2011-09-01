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

        // Get the result of the matching
        int nmatch = pm.getNumberOfMatchedVertices();
        int *mrows = new int[nmatch];
        int *mcols = new int[nmatch];
        int len;
        pm.getMatchedColumnsForRowsCopy(nmatch, len, mrows);
        pm.getMatchedRowsForColumnsCopy(nmatch, len, mcols);
        
        cout << endl << "Original Matrix:" << endl;
        std::cout<<*matrixPtr<<std::endl;

        // Create a new matrix with column permutation
        int max_entries = matrixPtr->MaxNumEntries();
        Epetra_CrsMatrix perm_matrix(Copy, matrixPtr->RowMap(), max_entries);
        int n =  matrixPtr->NumGlobalRows();

        double *values = new double[max_entries];
        int *indices = new int[max_entries];
        int num_entries;
        for (int i = 0; i < n ; i++)
        {
            // All in serial Comm so 0..n is fine
            matrixPtr->ExtractGlobalRowCopy(i, max_entries, num_entries, values,
                                        indices);
            for (int j = 0; j < num_entries; j++) indices[j]=mcols[indices[j]];
            perm_matrix.InsertGlobalValues(i, num_entries, values, indices);
        }
        perm_matrix.FillComplete();
        cout << endl << "After Column permutation:" << endl;
        cout << perm_matrix << endl;

        Epetra_CrsMatrix perm_matrix2(Copy, matrixPtr->RowMap(), max_entries);

        // Create a new matrix with row permutation
        for (int i = 0; i < n ; i++)
        {
            // All in serial Comm so 0..n is fine
            matrixPtr->ExtractGlobalRowCopy(i, max_entries, num_entries, values,
                                        indices);
            perm_matrix2.InsertGlobalValues(mrows[i], num_entries,
                                values, indices);
        }
        cout << endl << "After Row permutation:" << endl;
        perm_matrix2.FillComplete();
        cout << perm_matrix2 << endl;

        delete[] values;
        delete[] indices;
        delete[] mrows;
        delete[] mcols;

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
