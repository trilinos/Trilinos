#include "ifp_DenseMat.h"
/* Epetra_fmtflags typedef no longer used */
/*#include "Epetra_Object.h" // Bring in Epetra_fmtflags typedef */

ifp_DenseMat::ifp_DenseMat(const ifp_DenseMat& A)
{
    nrow = A.nrow; 
    ncol = A.ncol;
    register double *p = a = new double[nrow*ncol];
    register double *q = A.a;
    for (int i=0; i<nrow*ncol; i++)
        *p++ = *q++;
}

void ifp_DenseMat::Print(ostream& os) const
{
        // check not an implicit inverse
        assert (a != NULL);

/*        Epetra_fmtflags olda = os.setf(ios::right,ios::adjustfield);
        Epetra_fmtflags oldf = os.setf(ios::scientific,ios::floatfield);
        int oldp = os.precision(12); */

        const double *p = a;
        for (int j=0; j<numcol(); j++)
        for (int i=0; i<numrow(); i++)
               os << i+1 << "  " << j+1 << "  " << *p++ << endl;

/*        os.setf(olda,ios::adjustfield);
        os.setf(oldf,ios::floatfield);
        os.precision(oldp); */
}

// non-member functions

ostream& operator << (ostream& os, const ifp_DenseMat& mat)
{
        // should check not an implicit inverse

/*        Epetra_fmtflags olda = os.setf(ios::right,ios::adjustfield);
        Epetra_fmtflags oldf = os.setf(ios::scientific,ios::floatfield);
        int oldp = os.precision(12); */

        const double *a = &mat(0,0);
        for (int j=0; j<mat.numcol(); j++)
        for (int i=0; i<mat.numrow(); i++)
               os << i+1 << "  " << j+1 << "  " << *a++ << endl;

/*        os.setf(olda,ios::adjustfield);
        os.setf(oldf,ios::floatfield);
        os.precision(oldp); */
        return os;
}
