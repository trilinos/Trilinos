#include "ifp_LocalMat.h"
#include "ifp_BlockVec.h"
#include "ifp_BlockMat.h"
#include "ifp_biluk.h"
#include "ifp_ifpack.h"
#include "ifp_SparseUtil.h"
#include "ifp_DenseMat.h"
#include <iostream>
using namespace std;
/*int ifp_biluk::growth = 10; */

ifp_biluk::ifp_biluk()
{
    diag = (ifp_LocalMat **) NULL;
     al = NULL;
    jal = NULL;
    ial = NULL;
     au = NULL;
    jau = NULL;
    iau = NULL;
    NumEntries_ = 0;
    NumNonzeros_ = 0;
}

ifp_biluk::~ifp_biluk()
{
    int i;

    if (diag != NULL)
        for (i=0; i<Ap->numrow(); i++)
            delete diag[i];
    delete [] diag;

    if (al != NULL)
        for (i=0; i<ial[Ap->numrow()]; i++)
            delete al[i];
    delete [] al;

    if (au != NULL)
        for (i=0; i<iau[Ap->numrow()]; i++)
            delete au[i];
    delete [] au;

    delete [] jal;
    delete [] ial;
    delete [] jau;
    delete [] iau;
}

// if levfill < 0, then use existing factorization for the pattern

void ifp_biluk::setup(const ifp_BlockMat& A, int levfill)
{
    Ap = &A;

    if (levfill >= 0)
    {
	int nzl, nzu, ierr;

        delete jal;
        delete ial;
        delete jau;
        delete iau;
	int growthl = 2;
	int growthu = 2;
	ierr = 1;

	while (ierr !=0)
	  {
	    allocate_ilu(levfill, Ap->numrow(), &nzl, &nzu, 
			 &A.row_ptr(0), &A.col_ind(0), &ial, &jal, &iau, &jau, 
			 growthl, growthu);

	    ierr = symbolic_ilu(levfill, Ap->numrow(), &nzl, &nzu, 
				&A.row_ptr(0), &A.col_ind(0), ial, jal, iau, jau);
	    if (ierr !=0 )
	      {
		cout << "Doubling storage and trying again\n" << endl;
		delete [] ial;
		delete [] jal;
		delete [] iau;
		delete [] jau;
		if (ierr == -1) // Not enough for L
		  {
		    growthl *=2;
		    growthu *=2; // Assume current U size is too small
		  }
		else // ierr = -2
		  growthu *=2;
	      }
	  }
		    

        // allocate for values
  	delete [] al;
  	delete [] au;
        al = new ifp_LocalMatp[nzl];
        au = new ifp_LocalMatp[nzu];
    }

    int i, j, k, kk, id, idd;

    int *marker = new int[Ap->numrow()];
    for (i=0; i<Ap->numrow(); i++)
        marker[i] = -1; // negative flag

    // full length work array of matrices
    ifp_LocalMat **row = new ifp_LocalMatp[A.numcol()];  // array of pointers to matrices

    diag = new ifp_LocalMatp[A.numrow()];

    NumEntries_ = ial[A.numrow()] - ial[0] + iau[A.numrow()] - iau[0];
    NumNonzeros_ = 0;

    for (i=0; i<A.numrow(); i++)
    {
	int neqr = A.kvst_row(i+1) - A.kvst_row(i);

        // scatter data structure of L and U
	j = 0;
        for (k=ial[i]; k<ial[i+1]; k++)
	{
            marker[jal[k]] = j;
	    row[j] = A.val(0).CreateEmpty();
	    const int nvars = A.kvst_col(jal[k]+1)-A.kvst_col(jal[k]);
	    NumNonzeros_ += neqr*nvars;
	    row[j++]->SetToZero(neqr, nvars);
	}
        for (k=iau[i]; k<iau[i+1]; k++)
	{
            marker[jau[k]] = j;
	    row[j] = A.val(0).CreateEmpty();
	    const int nvars = A.kvst_col(jau[k]+1)-A.kvst_col(jau[k]);
	    NumNonzeros_ += neqr*nvars;
	    row[j++]->SetToZero(neqr, nvars);
	}
        
	//	cout << "============== Block Row "<< i << "==================" << endl;
	//	cout << "Number of rows in block row = "<< neqr << endl;
        // scatter row of A
        for (k=A.row_ptr(i); k<A.row_ptr(i+1); k++)
	  {
            row[marker[A.col_ind(k)]]->MatCopy(A.val(k));
	    /*
	    int nrow = ((ifp_DenseMat *)&(A.val(k)))->nrow;
	    int ncol = ((ifp_DenseMat *)&(A.val(k)))->ncol;
          cout << "Block Col "<< A.col_ind(k) << ",  Number of cols = " <<ncol << endl;
	    register double *qxx = ((ifp_DenseMat *)&(A.val(k)))->a;
	    for (int ixx=0; ixx<nrow; ixx++)
	      {
		for (int jxx=0; jxx<ncol; jxx++)
		  cout <<setiosflags(ios::scientific) <<*(qxx+jxx*nrow+ixx)<<",  ";
		cout << endl;
	      }
	    */
	  }
        ifp_LocalMat *mult = A.val(0).CreateEmpty();

        // eliminate the elements in L in order
        for (k=ial[i]; k<ial[i+1]; k++)
        {
            id = jal[k];

            // mult = row[id] / au[idiag[id]];
            row[marker[id]]->Mat_Mat_Mult(diag[id], mult, 1.0, 0.0);

            row[marker[id]]->MatCopy(*mult);

            for (kk=iau[id]+1; kk<iau[id+1]; kk++)
            {
                idd = jau[kk];
                if (marker[idd] >= 0) 
		{
                    // row[idd] = row[idd] - mult * au[kk];
                    mult->Mat_Mat_Mult(au[kk], row[marker[idd]], -1.0, 1.0);
		}
            }
        }

	delete mult;

        // gather resulting rows in L and U
        for (k=ial[i]; k<ial[i+1]; k++)
        {
            al[k] = row[marker[jal[k]]];
            marker[jal[k]] = -1;
        }
        for (k=iau[i]; k<iau[i+1]; k++)
        {
            au[k] = row[marker[jau[k]]];
            marker[jau[k]] = -1;
        }

        diag[i] = au[iau[i]]->CreateInv(local_precon);  // inverse of diagonal
    }

    delete [] row;
    delete [] marker;
}

void ifp_biluk::apply(int nr, int nc, const double *u, int ldu, 
    double *v, int ldv)
{
    applyl(nr, nc, u, ldu, v, ldv);
    applyr(nr, nc, v, ldv, v, ldv);
}

void ifp_biluk::applyl(int nr, int nc, const double *u, int ldu, 
    double *v, int ldv)
{
    int i, j;

    ifp_BlockVec  V(nr, nc, v, ldv, &Ap->kvst_col(0));
    ifp_BlockVec V2(nr, nc, v, ldv, &Ap->kvst_col(0));
    ifp_BlockVec  U(nr, nc, u, ldu, &Ap->kvst_col(0));
    V.VecCopy(U);

    // forward solve with lower triang factor (identity on diagonal assumed)

    for (i=0; i<Ap->numrow(); i++)
    {
        for (j=ial[i]; j<ial[i+1]; j++)
        {
            // V(i) = V(i) - al[j] * V(jal[j])
            al[j]->Mat_Vec_Mult(V(jal[j]), V2(i), -1.0, 1.0);
        }
    }

}

void ifp_biluk::applyr(int nr, int nc, const double *u, int ldu, 
    double *v, int ldv)
{
    int i, j;

    ifp_BlockVec  V(nr, nc, v, ldv, &Ap->kvst_col(0));
    ifp_BlockVec V2(nr, nc, v, ldv, &Ap->kvst_col(0));
    ifp_BlockVec  U(nr, nc, u, ldu, &Ap->kvst_col(0));
    V.VecCopy(U);

    // backward solve with upper triang factor

    for (i=Ap->numrow()-1; i>=0; i--)
    {
        for (j=iau[i]+1; j<iau[i+1]; j++)
        {
            // V(i) = V(i) - au[j] * V(jau[j])
            au[j]->Mat_Vec_Mult(V(jau[j]), V2(i), -1.0, 1.0);
        }
        diag[i]->Mat_Vec_Solve(V(i), V(i));
    }
}

void ifp_biluk::multiply(int nr, int nc, const double *u, int ldu, 
    double *v, int ldv)
{
    int i, j;

    ifp_BlockVec U(nr, nc, u, ldu, &Ap->kvst_col(0));
    ifp_BlockVec V(nr, nc, v, ldv, &Ap->kvst_col(0));
    ifp_BlockVec T(U);
    V.VecSetToZero();

    // multiply with upper triang factor

    for (i=0; i<Ap->numrow(); i++)
    {
        for (j=iau[i]; j<iau[i+1]; j++)
        {
            // V(i) = V(i) + au[j] * U(jau[j])
            au[j]->Mat_Vec_Mult(U(jau[j]), V(i), 1.0, 1.0);
        }
    }

    // multiply with lower triang factor (unit diagonal assumed)
    T.VecCopy(V);

    for (i=0; i<Ap->numrow(); i++)
    {
        for (j=ial[i]; j<ial[i+1]; j++)
        {
            // V(i) = V(i) + al[j] * temp(jal[j])
            al[j]->Mat_Vec_Mult(T(jal[j]), V(i), 1.0, 1.0);
        }
    }
}
double ifp_biluk::condest()
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
