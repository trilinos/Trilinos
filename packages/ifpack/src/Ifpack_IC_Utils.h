#ifndef IFPACK_IC_UTILS_H
#define IFPACK_IC_UTILS_H

typedef struct {
    double *val;  /* also known as A  */
    int    *col;  /* also known as JA; first column is column 0 */
    int    *ptr;  /* also known as IA; with ptr[0] = 0 */
} Ifpack_AIJMatrix;

extern "C" {
void quicksort (int *const pbase, double *const daux, int total_elems);
}

void Ifpack_AIJMatrix_dealloc(Ifpack_AIJMatrix *a);

void crout_ict(
    int n,
#ifdef IFPACK
    void * A,
    int maxentries,
    int (*getcol)( void * A, int col, int ** nentries, double * val, int * ind),
    int (*getdiag)( void *A, double * diag),
#else
    const Ifpack_AIJMatrix *AL,
    const double *Adiag,
#endif
    double droptol,
    int lfil,
    Ifpack_AIJMatrix *L,
    double **pdiag);


#endif
