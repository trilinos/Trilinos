#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <assert.h>
#include <string.h>
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Time.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

enum FillType {RHS_ONLY, MATRIX_ONLY, ALL}; 

class Fill { 
 public:
  void registerFillObjects( 
		  Epetra_Map &StandardMap, 
		  Epetra_Map &OverlapMap,
		  Epetra_Import & Importer,
		  Epetra_Comm &Comm);

  void fillMatrix(const Epetra_Vector *solnVector, Epetra_Vector *rhsVector, 
		  Epetra_RowMatrix *matrix);

 public:
  FillType flag;
  Epetra_Map *StandardMap; 
  Epetra_Map *OverlapMap;
  Epetra_Import *Importer;
  Epetra_Vector *soln;
  Epetra_CrsMatrix *A;
  Epetra_Vector *rhs;
  Epetra_Comm *Comm;
};



