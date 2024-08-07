// @HEADER
// *****************************************************************************
//  Zoltan Toolkit for Load-balancing, Partitioning, Ordering and Coloring
//
// Copyright 2012 NTESS and the Zoltan contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#include "phg.h"
#include "zz_const.h"
#include <limits.h>



#define SHOW_DISTMATRIX
/*
  #define SHOW_MINMAXV
  #define SHOW_MINMAXN
*/
#define SHOW_MINMAXP
    
    
char *Zoltan_PHG_uMe(PHGComm *hgc)
{
    static char msg[1024];

    sprintf(msg, "|%2d/%2d|: (%2d,%2d)/[%2d,%2d] ->", hgc->zz->Proc, hgc->zz->Num_Proc, hgc->myProc_x, hgc->myProc_y, hgc->nProc_x, hgc->nProc_y);
    return msg;
}

void Zoltan_PHG_uprintf(PHGComm *hgc, char *f_str,...)
{
va_list argp;

fflush(stdout);
fflush(stderr);
printf("%s", uMe(hgc)); 
va_start(argp, f_str);
vfprintf(stdout, f_str, argp);
va_end(argp);
fflush(stdout);
}

/*************************************************************************
* -------------------------- Error Exit ----------------------------------
**************************************************************************/
void Zoltan_PHG_errexit(char *f_str,...)
{
va_list argp;

fflush(stdout);
fflush(stderr);
fprintf(stderr, "\n****** Error:\n");
va_start(argp, f_str);
vfprintf(stderr, f_str, argp);
va_end(argp);

fprintf(stderr," ******\n");
fflush(stderr);
exit(1);
}


/*****************************************************************************/

void Zoltan_PHG_Find_Root(
  int val,
  int rank,
  MPI_Comm comm,
  int *bestval,
  int *bestrank
)
{
/* Based on local input value val, find the processor with the best val.
 * Return that processor and its value.
 * (Used when performing, say, local matching in each processor of a column 
 * and want to compute the best match in the column.)
 */
struct {
  int val;
  int rank;
} rootin, root;

    rootin.val = val;
    rootin.rank = rank;
    MPI_Allreduce(&rootin, &root, 1, MPI_2INT, MPI_MAXLOC, comm);

    *bestval = root.val;
    *bestrank = root.rank;
}


int Zoltan_PHG_LoadBalStat(ZZ *zz, HGraph *hg)
{
    char    *yo = "Zoltan_PHG_LoadBalStat";
    int     ierr=ZOLTAN_OK;
    PHGComm *comm = hg->comm;
    int     *v=NULL, *n=NULL, *p=NULL, x, y, i;
    int     minv=INT_MAX, maxv=-1, minn=INT_MAX, maxn=-1, minp=INT_MAX, maxp=-1;
    double  av=0.0, an=0.0, ap=0.0;

    if ((v = (int*) ZOLTAN_MALLOC(3 * comm->nProc * sizeof(int)))==NULL)
        MEMORY_ERROR;
    n = v + comm->nProc;
    p = n + comm->nProc;

    MPI_Gather(&hg->nVtx, 1, MPI_INT, v, 1, MPI_INT, 0, comm->Communicator);
    MPI_Gather(&hg->nEdge, 1, MPI_INT, n, 1, MPI_INT, 0, comm->Communicator);
    MPI_Gather(&hg->nPins, 1, MPI_INT, p, 1, MPI_INT, 0, comm->Communicator);

    for (i=0; i<comm->nProc; ++i) {
        minv = MIN(minv, v[i]);
        maxv = MAX(maxv, v[i]);
        av += v[i];
        minn = MIN(minn, n[i]);
        maxn = MAX(maxn, n[i]);
        an += n[i];
        minp = MIN(minp, p[i]);
        maxp = MAX(maxp, p[i]);
        ap += p[i];
    }

    av /= (double) comm->nProc;
    an /= (double) comm->nProc;
    ap /= (double) comm->nProc;
    
    if (!comm->myProc) {
#ifdef SHOW_DISTMATRIX        
        printf("Hypergraph distribution:\n     ");
        for (x=0; x<comm->nProc_x; ++x)
            printf("%-33d", x);
        printf("\n");
        for (y=0; y<comm->nProc_y; ++y) {
            printf("%3d: ", y);
            for (x=0; x<comm->nProc_x; ++x) {
                i = y* comm->nProc_x + x;
                printf("H(%7d, %7d, %9d)   ", v[i], n[i], p[i]);  
            }
            printf("\n");
            printf("     ");
            for (x=0; x<comm->nProc_x; ++x) {
                i = y* comm->nProc_x + x;
                printf("  ");
#ifdef SHOW_MINMAXV
                if (v[i]==minv)
                    printf("vvvvvvv  ");                
                else if (v[i]==maxv)
                    printf("^^^^^^^  ");                
                else
#endif
                    printf("         ");
#ifdef SHOW_MINMAXN                
                if (n[i]==minn)
                    printf("<<<<<<<  ");                
                else if (n[i]==maxn)
                    printf(">>>>>>>  ");                
                else
#endif
                    printf("         ");
#ifdef SHOW_MINMAXP                
                if (p[i]==minp)
                    printf("---------    ");                
                else if (p[i]==maxp)
                    printf("+++++++++    ");                
                else
#endif
                    printf("             ");
            }
            printf("\n");             
        }
#endif
        printf("Min:   (%7d, %7d, %9d)    Max: (%7d, %7d, %9d)\n", minv, minn, minp, maxv, maxn, maxp);
        printf("Imbal: (%7.2f, %7.2f, %9.2f)         (%7.2f, %7.2f, %9.2f)\n", 100.0*(av-minv)/av, 100.0*(an-minn)/an, 100.0*(ap-minp)/ap, 100.0*(maxv-av)/av, 100.0*(maxn-an)/an, 100.0*(maxp-ap)/ap);        
    }
 End:
    Zoltan_Multifree(__FILE__, __LINE__, 1, &v);
                         
    return ierr;
}

int Zoltan_PHG_isPrime(int n)
{
/* Naive program to test for primality. */
/* Returns accurate results for n <= maxValid. */
static const int maxValid = 250000;
static const int primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37,
                             41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97,
                             101, 103, 107, 109, 113, 127, 131, 137, 139, 
                             149, 151, 157, 163, 167, 173, 179, 181, 
                             191, 193, 197, 199, 211, 223, 227, 229,
                             233, 239, 241, 251, 257, 263, 269, 271, 277, 
                             281, 283, 293, 307, 311, 313, 317, 331, 337,
                             347, 349, 353, 359, 367, 373, 379, 383, 389,
                             397, 401, 409, 419, 421, 431, 433, 439,
                             443, 449, 457, 461, 463, 467, 479, 487, 491, 499};
static const int numprimes = sizeof(primes) / sizeof(int);
int i;
int rootn;
int isprime = 1;

 if (n == 1) return 0;

  rootn = sqrt((double)n)+1;
  for (i = 0; (i < numprimes) && (primes[i] < rootn); i++)
    if (!(n%primes[i])) {
      isprime = 0;
      break;
    }
  if (isprime && n>maxValid) {
    char str[128];  
    sprintf(str, "Warning: isPrime function may not be accurate for n(%i)>%d\n",
           n, maxValid);
    ZOLTAN_PRINT_WARN(-1, "Zoltan_PHG_isPrime", str);
  }
  
  return isprime;
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
