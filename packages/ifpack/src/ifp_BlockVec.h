/*@HEADER
// ************************************************************************
// 
//          Trilinos: An Object-Oriented Solver Framework
//              Copyright (2002) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//   
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//   
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
// 
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ************************************************************************
//@HEADER
*/

#ifndef _IFP_BLOCKVEC_H_
#define _IFP_BLOCKVEC_H_

#include <stdlib.h>
#ifndef NULL
#define NULL 0
#endif
#ifdef DEBUG
#include <assert.h>
#endif

class ifp_BlockVec
{
private:
    double *base;
    const int *partit;
    int owndata;

public:
    double *v;
    int     dim0;
    int     dim1;
    int     ld;
    int     size0;

    ifp_BlockVec& operator()(int i)
    {
#ifdef DEBUG
            assert(partit != NULL);
            assert(partit[i] < partit[i+1]);
            assert(partit[i+1] <= dim0);
#endif
	    v = base + partit[i];
	    size0 = partit[i+1] - partit[i];
	    return *this;
    }

    ifp_BlockVec(int nr, int nc, const double *a, int lda, 
                 const int *partitioning)
    {
        dim0 = nr;
        dim1 = nc;
        ld = lda;
	size0 = nr;
        base = (double *)a;
        v = (double *)a;
        partit = partitioning;
	owndata = 0;
    }

    ifp_BlockVec(const ifp_BlockVec& A);
    ifp_BlockVec(const ifp_BlockVec& A, int i);

   ~ifp_BlockVec()
    {
	if (owndata)
	    delete [] base;
	base = NULL;
    }

    void VecCopy(const ifp_BlockVec& A);
    void VecSetToZero();

    void BlockCopy(const ifp_BlockVec& A);
    void BlockSetToZero();
};

#endif // _IFP_BLOCKVEC_H_
