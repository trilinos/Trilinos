#ifndef _IFP_IFPACK_H_
#define _IFP_IFPACK_H_

#include "Ifpack_ConfigDefs.h"

#include <iostream>
#include <stdlib.h>
#include "ifp_arch.h"

#ifdef CRAY
extern "C" { void exit(int); }
#endif

#define ifp_error(msg,code) \
{ std::cerr << "IFPACK: " << msg << "    CODE " << code << std::endl; exit(1); }

#define MIN(x,y) ((x)<(y) ? (x) : (y))
#define MAX(x,y) ((x)>(y) ? (x) : (y))
#define ABS(x)   ((x)<0 ? (-(x)) : (x))
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TRUE
#define TRUE 1
#endif
#define SGN(x) ((x)<0.0 ? -1.0 : 1.0)

// miscellaneous prototypes for FORTRAN functions

extern "C"
{
void F77NAME(dgetri)(int *n, double *A, int *lda, int *ipiv,
    double *work, int *lwork, int *info);

void F77NAME(dgesvd) (char *, char *, int *, int *, const double *,
    int *, double *, double *, int *, double *, int *,
    double *, int *, int *);
}

#endif // _IFP_IFPACK_H_
