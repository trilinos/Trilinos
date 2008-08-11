/* Simple sparse matrix data structure. */

typedef struct {
    int i;     /* row index */
    int j;     /* column index */
    double val;  /* value */
} Matrix_entry;

typedef struct {
   int num_globalrows;  /* global # rows */
   int num_globalcols;  /* global # cols */
   int num_mynz;        /* local #nonzeros in matrix */
   Matrix_entry * entries; /* user must allocate space for the entries */
} Matrix;

int setup_zoltan(struct Zoltan_Struct *zz, Matrix A);
int run_zoltan(struct Zoltan_Struct *zz, Matrix A);

