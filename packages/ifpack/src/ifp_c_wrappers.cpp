#include "stdio.h"
#include "ifp_c_wrappers.h"
#include "ifp_BlockMat.h"
#include "ifp_brelax.h"
#include "ifp_biluk.h"
#include "az_aztec.h"
#include "ifp_ifpack.h"

#ifdef __cplusplus
using namespace std;
extern "C" {
#endif

// prototypes for this file
static void set_localprecon(ifp_GlobalPrecon *M, int local, 
  double lparam1, double lparam2);

// catch memory errors
//static void freeStoreException()
//{
//    cout << endl; // flush standard output
//    ifp_error("ifp_Fortran: free store exhausted", 0);
//}

// void ifp_initialize()
// {
//     void (*_new_handler)();
//     _new_handler = freeStoreException;
// }

void az2ifp_blockmatrix (void **bmat, AZ_MATRIX *Amat)
{

  double *val    = Amat->val;
  int *bindx  = Amat->bindx;
  int *indx = Amat->indx;
  int *bpntr = Amat->bpntr;
  int *rpntr = Amat->rpntr;
  int *cpntr = Amat->cpntr;
  int *data_org = Amat->data_org;
  int nrow = data_org[AZ_N_int_blk] + data_org[AZ_N_bord_blk];
  int ncol = nrow; /* Assume single processor for now; */
  int nnz = bpntr[nrow]; /* number of block nonzeros */
  int nnzs = indx[nnz]; /* number of scalar nonzeros */

  *bmat = (void *) new ifp_BlockMat(val, indx, bindx, rpntr, cpntr,
		   bpntr, nrow, ncol, nnz, nnzs);
}


void ifp_freeblockmatrix(void *bmat)
{
    delete (ifp_BlockMat *) bmat;
}

void ifp_preconditioner(
        void    **precon,
  const void    *bmat,
  const int     global,
  const double  gparam1,
  const double  gparam2,
  const int     local,
  const double  lparam1,
  const double  lparam2)
{
  ifp_GlobalPrecon *M;
  ifp_BlockMat *B = (ifp_BlockMat *) bmat;

  switch (global)
    {
    case 0:
      M = new ifp_None;
      break;
    case 1:
      M = new ifp_BJacobi;
      set_localprecon(M, (int)local, lparam1, lparam2);
      ((ifp_BJacobi *)M)->setup(*B);
      
      break;
    case 2:
      M = new ifp_BSOR;
      set_localprecon(M, (int)local, lparam1, lparam2);
      ((ifp_BSOR *)M)->setup(*B, (double)gparam1, (int)gparam2);
      break;
    case 3:
      M = new ifp_BSSOR;
      set_localprecon(M, (int)local, lparam1, lparam2);
      ((ifp_BSSOR *)M)->setup(*B, (double)gparam1, (int)gparam2);
      break;
    case 4:
      M = new ifp_biluk;
      set_localprecon(M, (int)local, lparam1, lparam2);
      ((ifp_biluk *)M)->setup(*B, (int)gparam1);
      
      break;
    default:
      ifp_error("ifp_c_wrapper: no such global preconditioner", 0);
    }
  *precon = (void *) M;
}

void ifp_freebiluk(void *precon)
{
    ifp_biluk * bilukPrec = (ifp_biluk *) precon;
    ifp_BlockMat * bilukMat = bilukPrec->A();
    delete bilukPrec;
    delete bilukMat;
}
void ifp_freepreconditioner(void *precon)
{
    delete (ifp_GlobalPrecon *) precon;
}

// static functions

#define INT(d) ((int)(d+0.5)) // round to integer

static void set_localprecon(ifp_GlobalPrecon *M, int local_, 
  double lparam1, double lparam2)
{
    LocalPreconName local = (LocalPreconName) local_;

    switch (local)
    {
    case LP_LU:
	M->localprecon(local);
	break;
    case LP_INVERSE:
	M->localprecon(local);
	break;
    case LP_SVD:
	M->localprecon(local, lparam1, lparam2);
	break;
    case LP_DIAG:
	M->localprecon(local);
	break;
    case LP_SOR:
	M->localprecon(local, lparam1, INT(lparam2));
	break;
    case LP_SSOR:
	M->localprecon(local, lparam1, INT(lparam2));
	break;
    default:
        ifp_error("ifp_Fortran: no such local preconditioner", 0);
    }
}

void ifp_matvec(
  void    *bmat,
  int     nr,
  int     nc,
  const double *u,
  int     ldu,
  double *v,
  int     ldv)
{
    ifp_BlockMat *B = (ifp_BlockMat *) bmat;
    B->mult(nr, nc, u, ldu, v, ldv);
}

void ifp_apply(
  void    *prec,
  int     nr,
  int     nc,
  const double *u,
  int     ldu,
  double *v,
  int     ldv)
{
    ifp_GlobalPrecon *M = (ifp_GlobalPrecon *) prec;
    M->apply(nr, nc, u, ldu, v, ldv);
}

void ifp_BJacobi_condest(void *M)
{
      double cond_number = ((ifp_BJacobi *)M)->condest();
      printf ("Condition number estimate = %.4e\n",cond_number);
}

void ifp_biluk_condest(void *M)
{
      double cond_number = ((ifp_biluk *)M)->condest();
      printf ("Condition number estimate = %.4e\n",cond_number);
}

void ifp_biluk_stats(void *M)
{
  ifp_biluk * Mp = (ifp_biluk * )M;
  ifp_BlockMat* Ap = Mp->A();
  int NumEntries_M = Mp->NumEntries();
  int NumNonzeros_M = Mp->NumNonzeros();
  int NumEntries_A = Ap->NumEntries();
  int NumNonzeros_A = Ap->NumNonzeros();
  double Entries_Ratio = ((double) NumEntries_M) / ((double) NumEntries_A);
  double Nonzeros_Ratio = ((double) NumNonzeros_M) / ((double) NumNonzeros_A);
  printf ("\n");
  printf ("********************************************************************\n");
  printf ("***** Number of Block Entries in User Matrix          = %d\n",NumEntries_A);
  printf ("***** Number of Nonzero Values in User Matrix         = %d\n",NumNonzeros_A);
  printf ("***** Number of Block Entries in IFPACK BILU Factors  = %d\n",NumEntries_M);
  printf ("***** Number of Nonzero Values in IFPACK BILU Factors = %d\n",NumNonzeros_M);
  printf ("***** Ratio of Block Entries in M to A                = %.4e\n",Entries_Ratio);
  printf ("***** Ratio of Nonzero Values in M to A               = %.4e\n",Nonzeros_Ratio);
  printf ("********************************************************************\n");
  printf ("\n");
}

#ifdef __cplusplus
}
#endif
