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

#include "ifp_BlockVec.h"

ifp_BlockVec::ifp_BlockVec(const ifp_BlockVec& A)
{
    dim0 = A.dim0;
    dim1 = A.dim1;
    ld = A.dim0;
    size0 = A.dim0;
    base = new double[dim0*dim1];
    v = base;
    partit = A.partit;
    owndata = 1;
    VecCopy(A);
}

ifp_BlockVec::ifp_BlockVec(const ifp_BlockVec& A, int index)
{
if (index >= 0)
{
    dim0 = A.partit[index+1] - A.partit[index];;
    dim1 = A.dim1;
    ld = dim0;
    size0 = dim0;
    base = new double[dim0*dim1];
    v = base;
    partit = NULL;
    owndata = 1;

    // copy the values

    int i, j;
    double *p;
    const double *q;

    for (i=0; i<dim1; i++)
    {
	p = v + i*ld;
        q = A.base + A.partit[index] + i*A.ld;
        for (j=0; j<size0; j++)
            *p++ = *q++;
    }
}
else
{
    dim0 = A.size0;
    dim1 = A.dim1;
    ld = dim0;
    size0 = dim0;
    base = new double[dim0*dim1];
    v = base;
    partit = NULL;
    owndata = 1;

    // copy the values

    int i, j;
    double *p;
    const double *q;

    for (i=0; i<dim1; i++)
    {
	p = v + i*ld;
        q = A.v + i*A.ld;
        for (j=0; j<size0; j++)
            *p++ = *q++;
    }
}
}

void ifp_BlockVec::VecCopy(const ifp_BlockVec& A)
{
#ifdef DEBUG
    assert(dim0 == A.dim0 && dim1 == A.dim1);
#endif
    int i, j;
    double *p;
    const double *q;

    for (i=0; i<dim1; i++)
    {
	p = v + i*ld;
        q = A.v + i*A.ld;
        for (j=0; j<dim0; j++)
            *p++ = *q++;
    }
}

void ifp_BlockVec::VecSetToZero()
{
    int i, j;
    double *p;

    for (i=0; i<dim1; i++)
    {
	p = v + i*ld;
        for (j=0; j<dim0; j++)
            *p++ = 0.0;
    }
}

void ifp_BlockVec::BlockCopy(const ifp_BlockVec& A)
{
#ifdef DEBUG
    assert(size0 == A.size0 && dim1 == A.dim1);
#endif
    int i, j;
    double *p;
    const double *q;

    for (i=0; i<dim1; i++)
    {
	p = v + i*ld;
        q = A.v + i*A.ld;
        for (j=0; j<size0; j++)
            *p++ = *q++;
    }
}

void ifp_BlockVec::BlockSetToZero()
{
    int i, j;
    double *p;

    for (i=0; i<dim1; i++)
    {
	p = v + i*ld;
        for (j=0; j<size0; j++)
            *p++ = 0.0;
    }
}
