#ifndef _ML_EPETRA_UTILS_H_
#define _ML_EPETRA_UTILS_H_

class Epetra_Comm;
class Epetra_BlockMap;
class Epetra_MultiVector;
class Epetra_RowMatrix;
#include "Epetra_LinearProblem.h"
#include "Epetra_Object.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Comm.h"
#include "Epetra_Import.h"

#ifdef ML_MPI
#ifndef EPETRA_MPI
#define EPETRA_MPI
#endif
#include "mpi.h"
#endif
#include "ml_include.h"

// External prototypes
extern "C" int Epetra_ML_matvec(void *data, int in, double *p, int out,
                 double *ap);

extern "C" int Epetra_ML_getrow(void *data, int N_requested_rows,
                 int requested_rows[], int allocated_space, int columns[],
                 double values[], int row_lengths[]);

extern "C" int Epetra_ML_comm_wrapper(double vec[], void *data);

extern "C" int EpetraMatrix2MLMatrix(ML *ml_handle, int level,
                                Epetra_RowMatrix * Amat);

extern "C" int Epetra2MLMatrix(Epetra_RowMatrix * A, ML_Operator *Result);

extern "C" Epetra_CrsMatrix *Epetra_MatrixMult(Epetra_RowMatrix *B, Epetra_RowMatrix *Bt);
extern "C" Epetra_CrsMatrix *Epetra_MatrixAdd(Epetra_RowMatrix *B, Epetra_RowMatrix *Bt, double scalar);
extern "C" int ML_Epetra_CRSinsert(ML_Operator *, int, int *, double *, int);


/* This Proto-type is in ml_rap.h. This prevents ml_rap.c and ml_matmat_mult.c     */
/* from including ml_epetra_utils.h which would require the C++ compiler for these */
/* files.             
extern int  ML_back_to_epetraCrs(ML_Operator *Mat1Mat2,  ML_Operator *Result, 
				 ML_Operator *Mat1,  ML_Operator *Mat2); 
*/



#endif /* _ML_EPETRA_UTILS_H_ */
