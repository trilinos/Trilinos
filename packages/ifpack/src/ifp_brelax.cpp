//@HEADER
// ***********************************************************************
// 
//       Ifpack: Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
//@HEADER

#include "ifp_LocalMat.h"
#include "ifp_BlockVec.h"
#include "ifp_BlockMat.h"
#include "ifp_brelax.h"
#include "ifp_ifpack.h"

void ifp_None::apply(int nr, int nc, const double *u, int ldu, 
    double *v, int ldv)
{
    int i, j;
    double *p;
    const double *q;

    for (i=0; i<nc; i++)
    {
        p = v + i*ldv;
        q = u + i*ldu;
        for (j=0; j<nr; j++)
            *p++ = *q++;
    }
}

ifp_BJacobi::ifp_BJacobi()
{
    Ap = (ifp_BlockMat *) NULL;
    diag = (ifp_LocalMat **) NULL;
}

ifp_BJacobi::~ifp_BJacobi()
{
    if (Ap != NULL)
        for (int i=0; i<Ap->numrow(); i++)
            delete diag[i];
    delete [] diag;
}

void ifp_BJacobi::setup(const ifp_BlockMat& A)
{
    int i, j;
    int got_diag;

    Ap = &A;
    diag = new ifp_LocalMatp[A.numrow()];

    // search for diagonal blocks
    for (i=0; i<A.numrow(); i++)
    {
        got_diag = FALSE;
        for (j=A.row_ptr(i); j<A.row_ptr(i+1); j++)
        {
            if (A.col_ind(j) == i)
            {
                diag[i] = A.val(j).CreateInv(local_precon);
                got_diag = TRUE;
            }
        }
        if (!got_diag)
            ifp_error("ifp_brelax: matrix does not have diagonal block", i);
    }
}

void ifp_BJacobi::apply(int nr, int nc, const double *u, int ldu,
    double *v, int ldv)
{
    ifp_BlockVec U(nr, nc, u, ldu, &Ap->kvst_row(0));
    ifp_BlockVec V(nr, nc, v, ldv, &Ap->kvst_row(0));

    for (int i=0; i<Ap->numrow(); i++)
    {
        diag[i]->Mat_Vec_Solve(U(i), V(i));
    }
}

double ifp_BJacobi::condest()
{
  // This routine computes a  bound of the infinity-norm condition number 
  // of the preconditioner.
  // It is equal to infinity-norm of M^-1 e, where e is the vector of all
  // ones.

  int i;

  int m = Ap->dimrow(); // Object is a matrix, so it has row/col dimensions.
  int n = Ap->dimcol();
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

// note:
// assumes each matrix row is stored such that lower triangular elements,
// diagonal elements, upper triangular elements are stored in order

ifp_BSOR_Base::ifp_BSOR_Base()
{
    Ap = (ifp_BlockMat *) NULL;
    diag = (ifp_LocalMat **) NULL;
    idiag = (int *) NULL;
}

ifp_BSOR_Base::~ifp_BSOR_Base()
{
    if (Ap != NULL)
        for (int i=0; i<Ap->numrow(); i++)
            delete diag[i];
    delete [] diag;
    delete [] idiag;
}

void ifp_BSOR_Base::setup(const ifp_BlockMat& A, double omega, int iterations)
{
    int i, j, nrow;
    int got_diag;

    omega_ = omega;
    iterations_ = iterations;
    Ap = &A;
    nrow = A.numrow();
    diag = new ifp_LocalMatp[nrow];
    idiag = new int[nrow];

    // search for diagonal blocks
    for (i=0; i<nrow; i++)
    {
        got_diag = FALSE;
        for (j=A.row_ptr(i); j<A.row_ptr(i+1); j++)
        {
            if (A.col_ind(j) == i)
            {
                diag[i] = A.val(j).CreateInv(local_precon);
                idiag[i] = j;
                got_diag = TRUE;
            }
        }
        if (!got_diag)
            ifp_error("ifp_brelax: matrix does not have diagonal block", i);
    }
}

// c = alpha * a + beta * b
// works on Blocks within a ifp_BlockVec

static void gaxpy(const double& alpha, const ifp_BlockVec& a,
    const double& beta, const ifp_BlockVec& b, ifp_BlockVec& c)
{
    double *ap, *bp, *cp;
    for (int j=0; j<a.dim1; j++)
    {
        ap = a.v + j * a.ld;
        bp = b.v + j * b.ld;
        cp = c.v + j * c.ld;

        for (int i=0; i<a.size0; i++) 
	    *cp++ = alpha * *ap++ + beta * *bp++;
            // c(i,j) = alpha*a(i,j) + beta*b(i,j);
    }
}

void ifp_BSOR::apply(int nr, int nc, const double *u, int ldu,
    double *v, int ldv)
{
    const int *ia = &Ap->row_ptr(0);
    const int *ja = &Ap->col_ind(0);
    int it, i, j;

#if 0
    ifp_BlockVec  V(nr, nc, v, ldv, &Ap->kvst_col(0));
    ifp_BlockVec V2(nr, nc, v, ldv, &Ap->kvst_col(0));
    ifp_BlockVec  U(nr, nc, u, ldu, &Ap->kvst_col(0));
    V.VecCopy(U);

    // Specialized code for first step.

    for (i=0; i<Ap->numrow(); i++)
    {
        for (j=ia[i]; j<idiag[i]; j++)
        {
            // V(i) = V(i) - omega_ a[j] * V(ja[j])
            Ap->val(j).Mat_Vec_Mult(V(ja[j]), V2(i), -omega_, 1.0);
        }
        diag[i]->Mat_Vec_Solve(V(i), V(i));
    }

    // After first step....
#endif
    ifp_BlockVec  V(nr, nc, v, ldv, &Ap->kvst_col(0));
    ifp_BlockVec  U(nr, nc, u, ldu, &Ap->kvst_col(0));
    V.VecSetToZero();

    for (it=1; it<=iterations_; it++)
    {
        for (i=0; i<Ap->numrow(); i++)
        {
            ifp_BlockVec temp(U,i);

            for (j=ia[i]; j<idiag[i]; j++)
            {
                // temp = temp - a[j] * v[ja[j]];
                Ap->val(j).Mat_Vec_Mult(V(ja[j]), temp, -1.0, 1.0);
            }
            for (j=idiag[i]+1; j<ia[i+1]; j++)
            {
                // temp = temp - a[j] * v[ja[j]];
                Ap->val(j).Mat_Vec_Mult(V(ja[j]), temp, -1.0, 1.0);
            }

            diag[i]->Mat_Vec_Solve(temp, temp);

            // v[i] = (1.0-omega_) * v[i] + omega_ * temp
            gaxpy(1.0-omega_, V(i), omega_, temp, V(i));
        }
    }
}

void ifp_BSSOR::apply(int nr, int nc, const double *u, int ldu,
    double *v, int ldv)
{
    const int *ia = &Ap->row_ptr(0);
    const int *ja = &Ap->col_ind(0);
    int it, i, j;

#if 0
    ifp_BlockVec  V(nr, nc, v, ldv, &Ap->kvst_col(0));
    ifp_BlockVec V2(nr, nc, v, ldv, &Ap->kvst_col(0));
    ifp_BlockVec  U(nr, nc, u, ldu, &Ap->kvst_col(0));
    V.VecCopy(U);

    // Specialized code for first step.

    // lower sweep
    for (i=0; i<Ap->numrow(); i++)
    {
        for (j=ia[i]; j<idiag[i]; j++)
        {
            // V(i) = V(i) - omega_ a[j] * V(ja[j])
            Ap->val(j).Mat_Vec_Mult(V(ja[j]), V2(i), -omega_, 1.0);
        }
        diag[i]->Mat_Vec_Solve(V(i), V(i));
    }

    // multiply by diagonal blocks
    for (i=0; i<Ap->numrow(); i++)
    {
        // V(i) = diag[i] * V(i)
        // this cannot be done in place
        ifp_BlockVec y(V,i); // make a copy of V(i)

        Ap->val(idiag[i]).Mat_Vec_Mult(y, V(i));
    }

    // upper sweep
    for (i=Ap->numrow()-1; i>=0; i--)
    {
        for (j=idiag[i]+1; j<ia[i+1]; j++)
        {
            // V(i) = V(i) - omega_ a[j] * V(ja[j])
            Ap->val(j).Mat_Vec_Mult(V(ja[j]), V2(i), -omega_, 1.0);
        }
        diag[i]->Mat_Vec_Solve(V(i), V(i));
    }

    // After first step....
#endif

    ifp_BlockVec  V(nr, nc, v, ldv, &Ap->kvst_col(0));
    ifp_BlockVec  U(nr, nc, u, ldu, &Ap->kvst_col(0));
    V.VecSetToZero();

    for (it=1; it<=iterations_; it++)
    {
        for (i=0; i<Ap->numrow(); i++)
        {
            ifp_BlockVec temp(U,i);

            for (j=ia[i]; j<idiag[i]; j++)
            {
                // temp = temp - a[j] * v[ja[j]];
                Ap->val(j).Mat_Vec_Mult(V(ja[j]), temp, -1.0, 1.0);
            }
            for (j=idiag[i]+1; j<ia[i+1]; j++)
            {
                // temp = temp - a[j] * v[ja[j]];
                Ap->val(j).Mat_Vec_Mult(V(ja[j]), temp, -1.0, 1.0);
            }

            diag[i]->Mat_Vec_Solve(temp, temp);

            // v[i] = (1.0-omega_) * v[i] + omega_ * temp
            gaxpy(1.0-omega_, V(i), omega_, temp, V(i));
        }

        for (i=Ap->numrow()-1; i>=0; i--)
        {
            ifp_BlockVec temp(U,i);

            for (j=ia[i]; j<idiag[i]; j++)
            {
                // temp = temp - a[j] * v[ja[j]];
                Ap->val(j).Mat_Vec_Mult(V(ja[j]), temp, -1.0, 1.0);
            }
            for (j=idiag[i]+1; j<ia[i+1]; j++)
            {
                // temp = temp - a[j] * v[ja[j]];
                Ap->val(j).Mat_Vec_Mult(V(ja[j]), temp, -1.0, 1.0);
            }

            diag[i]->Mat_Vec_Solve(temp, temp);

            // v[i] = (1.0-omega_) * v[i] + omega_ * temp
            gaxpy(1.0-omega_, V(i), omega_, temp, V(i));
        }
    }
}
