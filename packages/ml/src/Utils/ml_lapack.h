#ifndef AZTEC
#ifndef FSUB_TYPE
#  if defined(ncube)
#     define  FSUB_TYPE void
#  elif defined(paragon)
#     define  FSUB_TYPE void
#  elif defined(hp)
#     define  FSUB_TYPE void
#  else
#     define  FSUB_TYPE int
#  endif
#endif

extern FSUB_TYPE dgetrs_(char *, int *, int *, double *, int *, int *,
                          double *, int *, int *, unsigned int);

extern FSUB_TYPE  dgetrf_(int *, int *, double *, int *, int *, int *);

extern double ddot_(int *n1, double *v1, int *dum11, double *v2, int *dum21);

#else
#include "az_aztec.h"
#endif
extern FSUB_TYPE dgeqrf_(int *, int *, double *, int *,
        double *, double *, int *, int *);
extern FSUB_TYPE dorgqr_(int *, int *, int *, double *,
        int *, double *, double *, int *, int *);

