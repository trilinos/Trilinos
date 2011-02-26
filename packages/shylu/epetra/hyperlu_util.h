/** \file hyperlu_util.h

    \brief Utilities for HyperLU

    \author Siva Rajamanickam

*/

#include "Epetra_CrsMatrix.h" 

Epetra_CrsMatrix* balanceAndRedistribute(Epetra_CrsMatrix* A, 
                        Teuchos::ParameterList isoList);

void checkMaps(Epetra_CrsMatrix *A);

void findLocalColumns(Epetra_CrsMatrix *A, int *gvals, int &SNumGlobalCols);

void findBlockElems(int nrows, int *rows, int *gvals, int Lnr, int *LeftElems, 
        int Rnr, int *RightElems, string s1, string s2);
