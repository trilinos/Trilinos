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
