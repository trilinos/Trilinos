#include "ifp_BlockMat.h"
#include "ifp_BlockVec.h"
#include "ifp_DenseMat.h"
#include "ifp_ifpack.h"
#include "Epetra_Object.h" // Bring in Epetra_fmtflags typedef
/*#define DEBUG */

ifp_BlockMat::~ifp_BlockMat()
{
    if (a != NULL)
    {
      //delete [] a[0]->Data();       // values
            for (int i=0; i<nnz; i++)
                a[i]->Data() = NULL;      // so DenseMat destructor safe
            delete [] (ifp_DenseMat *) a[0];  // must be array of ifp_DenseMat
    }
    delete [] a;
    //delete [] ja;
    //delete [] ia;
    //delete kvstr;
    //delete kvstc;
}

void ifp_BlockMat::mult(int nr, int nc, const double *B, int ldu, 
    double *C, int ldv) const
{
    ifp_BlockVec Bb(nr, nc, B, ldu, kvstc);
    ifp_BlockVec Cb(dimrow(), nc, C, ldv, kvstr);
    Cb.VecSetToZero();

    for (int i=0; i<nrow; i++)
    {
        for (int j=ia[i]; j<ia[i+1]; j++)
	{
             a[j]->Mat_Vec_Mult(Bb(ja[j]), Cb(i), 1.0, 1.0);
	}
    }
}

void ifp_BlockMat::trans_mult(int nr, int nc, const double *B, int ldu,
    double *C, int ldv) const
{
    ifp_BlockVec Bb(nr, nc, B, ldu, kvstc);
    ifp_BlockVec Cb(dimcol(), nc, C, ldv, kvstr);
    Cb.VecSetToZero();

    for (int i=0; i<nrow; i++)
    {
        for (int j=ia[i]; j<ia[i+1]; j++)
	{
             a[j]->Mat_Trans_Vec_Mult(Bb(i), Cb(ja[j]), 1.0, 1.0);
	}
    }
}

// Converts Aztec Matrix to IFPACK matrix

ifp_BlockMat::ifp_BlockMat(double *val, int *indx, int * bindx, int *rpntr, int *cpntr,
		   int *bpntr, int nrow_, int ncol_, int nnz_, int nnzs_)
{
    int i, j, k;
    double *pp;

    a = (ifp_LocalMat **) NULL;
    ja = bindx;
    ia = bpntr;
    kvstr = rpntr;
    kvstc = cpntr;
    nrow = nrow_;
    ncol = ncol_;
    nnz = nnz_;
    nnzs = nnzs_;

#ifdef DEBUG
     cerr << nrow << " block rows.\n";
     cerr << ncol << " block columns.\n";

     cerr << nnz << " block nonzeros.\n";
     cerr << nnzs << " stored scalar entries.\n";
#endif

    a  = new ifp_LocalMatp[nnz];   // array of pointers to blocks

    ifp_DenseMat *p = new ifp_DenseMat[nnz]; // array of blocks

    for (k=0; k<nnz; k++)
        a[k] = &p[k];

    pp = val;  // array of values

    // fill entries and ja -----------------------------------------------------

    k = 0;

    // loop on block rows
    for (i=0; i<nrow; i++)
    {
        int neqr = kvstr[i+1]-kvstr[i];

	for (j=ia[i]; j<ia[i+1]; j++)
        {
	  int neqc = kvstc[bindx[j]+1] - kvstc[bindx[j]];
	  p[k].set(neqr,neqc,pp);
	  /*	  printf("In az2bp i = %d  j = %d\n",i,j);
		  cout << p[k]; */
	  pp = pp + neqr * neqc;
	  k++;
        }
    }

}

// Non-member functions

ostream& operator << (ostream& os, const ifp_BlockMat& mat)
{
        int M = mat.numrow();
        int N = mat.numcol();
        int rowp1, colp1;
        int flag = 0;
        Epetra_fmtflags olda = os.setf(ios::right,ios::adjustfield);
        Epetra_fmtflags oldf = os.setf(ios::scientific,ios::floatfield);
        int oldp = os.precision(12);

        for (int i = 0; i < M ; i++)
           for (int j=mat.row_ptr(i);j<mat.row_ptr(i+1);j++)
           {   
              rowp1 =  i + 1;
              colp1 =  mat.col_ind(j) + 1;
              if ( rowp1 == M && colp1 == N ) flag = 1;
              os.width(14);
              os <<  rowp1 ; os << "    " ;
              os.width(14);
              os <<  colp1;
#if 0
              os.width(20);
              os <<  *(mat.val(j).Data());
#endif
              os << endl;
           }

#if 0
        if (flag == 0)
        {
           os.width(14);
           os <<  M ; os << "    " ;
           os.width(14);
           os <<  N ; os << "    " ;
           os.width(20);
           os <<  mat(M-1,N-1) << "\n";
        }
#endif

        os.setf(olda,ios::adjustfield);
        os.setf(oldf,ios::floatfield);
        os.precision(oldp);

        return os;
}
