#include "ifp_Matrix.h"
#include "ifp_Precon.h"
#include "ifp_ifpack.h"

double ifp_Precon::condest()
{
  // This routine computes a  bound of the infinity-norm condition number 
  // of the preconditioner.
  // It is equal to infinity-norm of M^-1 e, where e is the vector of all
  // ones.

  int i;

  int m = dimrow(); // Object is a matrix, so it has row/col dimensions.
  int n = dimcol();
  double *u = new double[n];
  double *v = new double[m];

    for (i=0; i<n; i++) u[i] = 1.0;

    apply(n, 1, u, n, v, n);

    double cond_number = 0.0;
    for (i=0; i<m; i++)
        if (ABS(v[i]) > cond_number)
           cond_number = ABS(v[i]);

    delete [] u;
    delete [] v;

    return(cond_number);
}
