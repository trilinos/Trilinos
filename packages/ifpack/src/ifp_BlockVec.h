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
