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
