#include "ifp_ifpack.h"
#include "ifp_SparseUtil.h"
#include "stdio.h" // kludge
using namespace std;
/* #define DEBUG */
// shell sort
// stable, so it is fast if already sorted
// sorts x[0:n-1] in place, ascending order.

void shell_sort(
  const int n,
  int x[])
{
    int m, max, j, k, itemp;
    
    m = n/2;

    while (m > 0) {
        max = n - m;
        for (j=0; j<max; j++)
        {
            for (k=j; k>=0; k-=m)
            {
                if (x[k+m] >= x[k])
                    break;
                itemp = x[k+m];
                x[k+m] = x[k];
                x[k] = itemp;
            }
        }   
        m = m/2;
    }
}

// allocate space for integer data for level ILU

void allocate_ilu(
  const int levfill,                 // level of fill
  const int n,                       // order of matrix
  int *nzl, int *nzu,                // space allocated
  const int ia[], const int ja[],    // input
  int *ial[], int *jal[],            // output lower factor structure
  int *iau[], int *jau[],            // output upper factor structure
  int growthl, int growthu)              // storage parameters
{
    int i, j, nzla, nzua;
    int maxl, maxu;

    maxu = n*(n+1)/2;
    maxl = n*n - maxu;

    // count number of entries in lower and upper triangular parts

    nzla = 0;
    nzua = 0;
    for (i=0; i<n; i++)
    {
        for (j=ia[i]; j<ia[i+1]; j++)
        {
            if (ja[j] < i)
               nzla++;
            else
               nzua++;
        }
    }

    *ial = new int[n+1];
    *iau = new int[n+1];

    if (levfill == 0) // ILU(0)
    {
        *nzl = nzla;
        *nzu = nzua + n;
    }
    else if (levfill == n) // full factorization
    {
        *nzl = maxl;
        *nzu = maxu;
    }
    else
    {
        *nzl = MIN((levfill+growthl)*nzla, maxl);
        *nzu = MIN((levfill+growthu)*nzua + n, maxu);
    }

    *jal = new int[*nzl];
    *jau = new int[*nzu];

    if (levfill != 0)
    {
      // cerr << "nnz in unfactored mat: " << nzla + nzua << endl;
      // cerr << "nnz allocated for ILU: " << *nzl + *nzu << endl;
    }
}

// symbolic level ILU
// fortran-style data structures and algorithms
// factors into separate upper and lower parts
// assumes rows are sorted  *** symbolic ILU
// assumes no zero rows
// Returns an integer:
//  0 - No error
// -1 - Not enough memory for L
// -2 - Not enough memory for U

int symbolic_ilu(
  const int levfill,                 // level of fill
  const int n,                       // order of matrix
  int *nzl,                          // input-output
  int *nzu,                          // input-output
  const int ia[], const int ja[],    // input
  int ial[], int jal[],              // output lower factor structure
  int iau[], int jau[])              // output upper factor structure
{
#if 0
    // copy if levfill is 0
    if (levfill == 0)
    {
        int kl = 0;
        int ku = 0;
        ial[0] = 0;
        iau[0] = 0;
        for (int i=0; i<n; i++)
        {
            for (int j=ia[i]; j<ia[i+1]; j++)
            {
                if (ja[j] < i)
                   jal[kl++] = ja[j];
                else
                   jau[ku++] = ja[j];
            }
            ial[i+1] = kl;
            iau[i+1] = ku;
            shell_sort(ial[i+1]-ial[i], &jal[ial[i]]);
            shell_sort(iau[i+1]-iau[i], &jau[iau[i]]);
        }
        assert (kl <= *nzl); 
        assert (ku <= *nzu);
        *nzl = kl;
        *nzu = ku;
    }
    else
#endif
    {

    int *lnklst = new int[n];
    int *curlev = new int[n];
    int *levels = new int[*nzu];
    int *iwork = new int[n];

    int knzl = 0;
    int knzu = 0;

    ial[0] = 0;
    iau[0] = 0;

    for (int i=0; i<n; i++)
    {
        int first, next, j;

        // copy column indices of row into workspace and sort them

        int len = ia[i+1] - ia[i];
	if (len == 0 ) cerr << " Row " << i << " is zero " << endl;
        next = 0;
        for (j=ia[i]; j<ia[i+1]; j++)
            iwork[next++] = ja[j];
        shell_sort(len, iwork);

        // construct implied linked list for row

        first = iwork[0];
        curlev[first] = 0;

        for (j=0; j<=len-2; j++)
        {
            lnklst[iwork[j]] = iwork[j+1];
            curlev[iwork[j]] = 0;
        }

        lnklst[iwork[len-1]] = n;
        curlev[iwork[len-1]] = 0;

        // merge with rows in U

        next = first;
        while (next < i)
        {
            int oldlst = next;
            int nxtlst = lnklst[next];
            int row = next;
            int ii;

            // scan row

            for (ii=iau[row]+1; ii<iau[row+1]; /*nop*/)
            {
                if (jau[ii] < nxtlst)
                {
                    // new fill-in
                    int newlev = curlev[row] + levels[ii] + 1;
                    if (newlev <= levfill)
                    {
                        lnklst[oldlst]  = jau[ii];
                        lnklst[jau[ii]] = nxtlst;
                        oldlst = jau[ii];
                        curlev[jau[ii]] = newlev;
                    }
                    ii++;
                }
                else if (jau[ii] == nxtlst)
                {
                    oldlst = nxtlst;
                    nxtlst = lnklst[oldlst];
                    int newlev = curlev[row] + levels[ii] + 1;
                    curlev[jau[ii]] = MIN(curlev[jau[ii]], newlev);
                    ii++;
                }
                else // (jau[ii] > nxtlst)
                {
                    oldlst = nxtlst;
                    nxtlst = lnklst[oldlst];
                }
            }
            next = lnklst[next];
        }
        
        // gather the pattern into L and U

        next = first;
        while (next < i)
        {
           if (knzl >= *nzl)
	     {
	       cerr << "For row "<< i << endl;
	       cerr << " nzl = "<< *nzl <<" knzl = "<< knzl << endl;
	       printf("Not enough space allocated for lower symbolic factor\n");
	       return(-1);
	     }

            jal[knzl++] = next;
            next = lnklst[next];
        }
        ial[i+1] = knzl;

        if (next != i)
        {
            cerr << i << "  U has zero on diag, forcing nonzero" << endl;
            if (knzu >= *nzu)
	     {
	       cerr << " nzl = "<< *nzu <<" knzl = "<< knzu << endl;
	       printf("Not enough space allocated for symbolic factor\n");
	       return(-2);
	     }
            levels[knzu] = 2*n; // infinity
            jau[knzu++] = i;
        }

        while (next < n)
        {
            if (knzu >= *nzu)
	     {
	       cerr << "For row "<< i << endl;
	       cerr << " nzu = "<< *nzu <<" knzu = "<< knzu << endl;
	       printf("Not enough space allocated for upper symbolic factor\n");
	       return(-2);
    
	     }
            levels[knzu] = curlev[next];
            jau[knzu++] = next;
            next = lnklst[next];
        }
        iau[i+1] = knzu;
    }

    delete [] lnklst;
    delete [] curlev;
    delete [] levels;
    delete [] iwork;

    *nzl = knzl;
    *nzu = knzu;

#ifdef DEBUG
    cerr << "Actual nnz for ILU: " << *nzl + *nzu << endl;
#endif
    }
    return(0);
}
