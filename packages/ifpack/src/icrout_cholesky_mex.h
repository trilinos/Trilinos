
#define IFPACK
/*
 * Data structure for sparse matrices is CSR, 0-based indexing.
 */
typedef struct {
    double *val;  /* also known as A  */
    int    *col;  /* also known as JA; first column is column 0 */
    int    *ptr;  /* also known as IA; with ptr[0] = 0 */
} Matrix;

void Matrix_dealloc(Matrix *a);
void crout_ict(
    int n,
#ifdef IFPACK
    void * A,
    int maxentries;
    int (*getcol)( void * A, int col, int ** nentries, double * val, int * ind),
    int (*getdiag)( void *A, double * diag),
#else
    const Matrix *AL,
    const double *Adiag,
#endif
    double droptol,
    int lfil,
    Matrix *L,
    double **pdiag);
