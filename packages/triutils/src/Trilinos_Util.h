class Epetra_Comm;
class Epetra_Vector;
class Epetra_Map;
class Epetra_BlockMap;
class Epetra_CrsMatrix;
class Epetra_VbrMatrix;
class Epetra_MultiVector;

#define Trilinos_Util_max(x,y) (( x > y ) ? x : y)     /* max function  */
#define Trilinos_Util_min(x,y) (( x < y ) ? x : y)     /* min function */

#ifdef HAVE_CONFIG_H
#include "Triutils_config.h"

#ifdef HAVE_DEBUG
#ifndef DEBUG
#define DEBUG
#endif
#endif

#endif /*HAVE_CONFIG_H*/

void Trilinos_Util_read_hb(char *data_file, int MyPID,
	      int *N_global, int *n_nonzeros, 
	      double **val, int **bindx,
	      double **x, double **b, double **xexact);

void Trilinos_Util_read_coo(char *data_file, int MyPID,
	      int *N_global, int *n_nonzeros,
	      double **val, int **bindx,
	      double **x, double **b, double **xexact);

void Trilinos_Util_ReadHb2Epetra(char *data_file,
				 const Epetra_Comm  &comm, 
				 Epetra_Map *& map, 
				 Epetra_CrsMatrix *& A, 
				 Epetra_Vector *& x, 
				 Epetra_Vector *& b,
				 Epetra_Vector *&xexact);

void Trilinos_Util_ReadHpc2Epetra(char *data_file,
				 const Epetra_Comm  &comm, 
				 Epetra_Map *& map, 
				 Epetra_CrsMatrix *& A, 
				 Epetra_Vector *& x, 
				 Epetra_Vector *& b,
				 Epetra_Vector *&xexact);

void Trilinos_Util_ReadHb2EpetraVbr(char *data_file, char * partitioning,
				 const Epetra_Comm  &comm, 
				 Epetra_BlockMap *& map, 
				 Epetra_VbrMatrix *& A, 
				 Epetra_Vector *& x, 
				 Epetra_Vector *& b,
				 Epetra_Vector *&xexact);

void Trilinos_Util_distrib_msr_matrix(const Epetra_Comm & Comm,
	      int *N_global, int *n_nonzeros, 
           int *N_update, int **update, 
	      double **val, int **bindx,
	      double **x, double **b, double **xexact);

void Trilinos_Util_distrib_vbr_matrix(const Epetra_Comm & Comm,
	      int *N_global, int *N_blk_global, 
           int *n_nonzeros,  int *n_blk_nonzeros,
           int *N_update, int **update, 
	      double **val, int **indx, int **rpntr, int **cpntr,
           int **bpntr, int **bindx,
	      double **x, double **b, double **xexact);

void Trilinos_Util_create_vbr(const Epetra_Comm & Comm, char *part_file,
                int *N_global, int *N_blk_global, 
                int *n_nonzeros, int *n_blk_nonzeros,
                int *N_update, int **update,
		      int *bindx_msr, double *val_msr,
		      double **val, int **indx, int **rpntr, int **cpntr,
		      int **bpntr, int **bindx);

double Trilinos_Util_smsrres (int m, int n, 
	      double *val, int *indx, 
	      double *xlocal, double *x, double *b);

double Trilinos_Util_scscres (int isym, int m, int n, 
	      double *val, int *indx, int *pntr,
	      double *x, double *b);

void  Trilinos_Util_scscmv (int isym, int m, int n, 
	      double *val, int *indx, int *pntr,
	      double *x, double *b);

double Trilinos_Util_svbrres (int m, int n, int m_blk,
		double *val, int *indx, int *bindx, int *rpntr,
		int *cpntr, int *bpntrb, int *bpntre,
		double *x, double *b);

void Trilinos_Util_msr2vbr(double val[], int indx[], int rnptr[], int
                cnptr[], int bnptr[],
                int bindx[], int msr_bindx[], double msr_val[],
                int total_blk_rows, int total_blk_cols, int blk_space,
                int nz_space, int blk_type);

int Trilinos_Util_find_block_col(int cnptr[], int column, int
                                 max_blocks, int blk_size);

int Trilinos_Util_find_block_in_row(int bindx[], int bnptr[], int
                                    blk_row, int blk_col,
                                    int indx[], int no_elements, double val[],
                                    int blk_space, int nz_space);


void Trilinos_Util_add_new_ele(int cnptr[], int col, int blk_row, int
                               bindx[], int bnptr[],
                    int indx[], double val[], int row, double new_ele,
                    int maxcols, int blk_space, int nz_space, int blk_type);

int Trilinos_Util_find_closest_not_larger(int key, int list[], int length);

void Trilinos_Util_convert_values_to_ptrs(int array[], int length, int start);

int Trilinos_Util_csrcsc(int n, int n2, int job, int ipos, double * a,
           int *ja, int *ia, double *ao, int *jao, int *iao);

int Trilinos_Util_csrmsr( int n, double *a, int *ja, int *ia, double *ao,
             int *jao, double *wk, int *iwk);

int Trilinos_Util_ssrcsr( int job, int value2, int nrow, double *a,
                          int *ja, int *ia, int nzmax,
                          double *ao, int *jao, int *iao, int *indu, 
                          int *iwk);

int Trilinos_Util_coocsr( int nrow, int nnz, double *a, int *ir, int *jc,
                          double *ao, int *jao, int *iao);


struct SPBLASMAT_STRUCT {
  int n;
  double *val;
  int *indx;
  int *bindx;
  int *rpntr;
  int *cpntr;
  int *bpntrb;
  int *bpntre;
  int buffersize;
  int bufferstride;
  double *buffer;
  int *ncolvec;
  double nops_per_rhs;
  int minblocksize;
  int maxblocksize;
};
typedef struct SPBLASMAT_STRUCT SPBLASMAT;

#define MAXNRHS 1

void  Trilinos_Util_duscr_vbr(int n, double *val, int *indx, int *bindx, 
		      int *rpntr, int *cpntr, int *bpntrb, int *bpntre, 
		      SPBLASMAT *A);
void Trilinos_Util_dusmm(int m, int nrhs, int k, double alpha, SPBLASMAT *A,
		 double *x, int xstride, double beta, double *b, int bstride);
void  Trilinos_Util_dusds_vbr( SPBLASMAT *A);

void Trilinos_Util_GenerateCrsProblem(int nx, int ny, int npoints, int * xoff, int * yoff,
																			const Epetra_Comm  &comm, 
																			Epetra_Map *& map, 
																			Epetra_CrsMatrix *& A, 
																			Epetra_Vector *& x, 
																			Epetra_Vector *& b,
																			Epetra_Vector *&xexact);

void Trilinos_Util_GenerateCrsProblem(int nx, int ny, int npoints, int * xoff, int * yoff, int nrhs,
																			const Epetra_Comm  &comm, 
																			Epetra_Map *& map, 
																			Epetra_CrsMatrix *& A, 
																			Epetra_MultiVector *& x, 
																			Epetra_MultiVector *& b,
																			Epetra_MultiVector *&xexact);

