
/*
 * File name:	pzlangs.c
 * History:     Modified from lapack routine ZLANGE
 */
#include <math.h>
#include "superlu_zdefs.h"

double pzlangs(char *norm, SuperMatrix *A, gridinfo_t *grid)
{
/* 
    Purpose   
    =======   

    PZLANGS returns the value of the one norm, or the Frobenius norm, or 
    the infinity norm, or the element of largest absolute value of a 
    real matrix A.   

    Description   
    ===========   

    PZLANGE returns the value   

       PZLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'   
                 (   
                 ( norm1(A),         NORM = '1', 'O' or 'o'   
                 (   
                 ( normI(A),         NORM = 'I' or 'i'   
                 (   
                 ( normF(A),         NORM = 'F', 'f', 'E' or 'e'   

    where  norm1  denotes the  one norm of a matrix (maximum column sum), 
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and 
    normF  denotes the  Frobenius norm of a matrix (square root of sum of 
    squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.   

    Arguments   
    =========   

    NORM    (input) CHARACTER*1   
            Specifies the value to be returned in DLANGE as described above.   
    A       (input) SuperMatrix*
            The M by N sparse matrix A. 
    GRID    (input) gridinof_t*
            The 2D process mesh.
   ===================================================================== 
*/
    
    /* Local variables */
    NRformat_loc *Astore;
    int_t    m_loc;
    doublecomplex   *Aval;
    int_t    i, j, irow, jcol;
    double   value=0., sum;
    double   *rwork;
    double   tempvalue;
    double   *temprwork;

    Astore = (NRformat_loc *) A->Store;
    m_loc = Astore->m_loc;
    Aval   = (doublecomplex *) Astore->nzval;
    
    if ( SUPERLU_MIN(A->nrow, A->ncol) == 0) {
	value = 0.;
    } else if (lsame_(norm, "M")) {
	/* Find max(abs(A(i,j))). */
	value = 0.;
	for (i = 0; i < m_loc; ++i) {
	    for (j = Astore->rowptr[i]; j < Astore->rowptr[i+1]; ++j)
		value = SUPERLU_MAX( value, z_abs(&Aval[j]) );
	}

	MPI_Allreduce(&value, &tempvalue, 1, MPI_DOUBLE, MPI_MAX, grid->comm);
	value = tempvalue;

    } else if (lsame_(norm, "O") || *(unsigned char *)norm == '1') {
	/* Find norm1(A). */
	value = 0.;
#if 0
	for (j = 0; j < A->ncol; ++j) {
	    sum = 0.;
	    for (i = Astore->colptr[j]; i < Astore->colptr[j+1]; i++) 
		sum += fabs(Aval[i]);
	    value = SUPERLU_MAX(value,sum);
	}
#else /* XSL ==> */
	if ( !(rwork = (double *) doubleCalloc_dist(A->ncol)) )
	    ABORT("doubleCalloc_dist fails for rwork.");
	for (i = 0; i < m_loc; ++i) {
	    for (j = Astore->rowptr[i]; j < Astore->rowptr[i+1]; ++j) {
	        jcol = Astore->colind[j];
		rwork[jcol] += z_abs(&Aval[j]);
	    }
	}

	if ( !(temprwork = (double *) doubleCalloc_dist(A->ncol)) )
	    ABORT("doubleCalloc_dist fails for temprwork.");
	MPI_Allreduce(rwork, temprwork, A->ncol, MPI_DOUBLE, MPI_SUM, grid->comm);
	value = 0.;
	for (j = 0; j < A->ncol; ++j) {
	    value = SUPERLU_MAX(value, temprwork[j]);
	}
	SUPERLU_FREE (temprwork);
	SUPERLU_FREE (rwork);
#endif	
    } else if (lsame_(norm, "I")) {
	/* Find normI(A). */
	value = 0.;
	sum = 0.;
	for (i = 0; i < m_loc; ++i) {
	    for (j = Astore->rowptr[i]; j < Astore->rowptr[i+1]; ++j)
	        sum += z_abs(&Aval[j]);
	    value = SUPERLU_MAX(value, sum);
	}
	MPI_Allreduce(&value, &tempvalue, 1, MPI_DOUBLE, MPI_MAX, grid->comm);
	value = tempvalue;

    } else if (lsame_(norm, "F") || lsame_(norm, "E")) {
	/* Find normF(A). */
	ABORT("Not implemented.");
    } else {
	ABORT("Illegal norm specified.");
    }
    
    return (value);

} /* pzlangs */
