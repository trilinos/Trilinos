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

#define MAXNRHS 16

#define name2(a,b) a ## b

#ifndef HAVE_NO_FORTRAN_UNDERSCORE
#define F77NAME(x) name2(x,_)
#else
#define F77NAME(x) x
#endif

#ifdef ultra
#define _TIMER second_linux
#endif
#ifdef alpha
#define _TIMER second_alpha
#endif
#ifdef rs6000
#define _TIMER F77NAME(second_rs6000)
#endif
#ifdef pclinux
#define _TIMER second_linux
#endif
#ifdef pgilinux
#define _TIMER second_linux
#endif
#ifdef irix64
#define _TIMER F77NAME(second_rs6000)
#endif

double F77NAME(dnrm2)(int *n, double *x, int *incx);
void   F77NAME(d27ptgen)
     (int *nnzmx, int *nx, int *ny, int *nz, int *nb, 
      double *factor, int *izero, char *title, int *n, 
      double *val, int *bpntr, int *bindx, int *nnz, int *ierr);
void  cblas_duscr_vbr(int n, double *val, int *indx, int *bindx, 
		      int *rpntr, int *cpntr, int *bpntrb, int *bpntre, 
		      SPBLASMAT *A);
void cblas_dusmm_ref(int m, int nrhs, int k, double alpha, SPBLASMAT *A,
		 double *x, int xstride, double beta, double *b, int bstride);

void cblas_dusmm(int m, int nrhs, int k, double alpha, SPBLASMAT *A,
		 double *x, int xstride, double beta, double *b, int bstride);

void F77NAME(daxpy) (int *n, double *alpha, double *x, 
		     int *incx, double *y, int *incy);
void F77NAME(dgemm) ( char *TRANSA, char *TRANSB, int *M, int *N, int *K,
                      double *ALPHA, double *A, int *LDA, double *B, int *LDB,
                      double *BETA, double *C, int *LDC);

void _TIMER(double *seconds);
#define MAX(a,b)  ( ( (a) > (b) ) ? (a) : (b) )
#define MIN(a,b)  ( ( (a) < (b) ) ? (a) : (b) )
