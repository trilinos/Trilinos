/*
** $Id$
**
** Miscellaneous functions to support writing simple Zoltan examples.
*/
#include "exzoltan.h"
#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif


void exShowError(int val, char *s, int me)
{
  if (s)
    {
    printf("%s ",s);
    }
  switch (val)
    {
    case ZOLTAN_OK:
      printf("%d: SUCCESSFUL\n", me);
      break;
    case ZOLTAN_WARN:
      printf("%d: WARNING\n", me);
      break;
    case ZOLTAN_FATAL:
      printf("%d: FATAL ERROR\n", me);
      break;
    case ZOLTAN_MEMERR:
      printf("%d: MEMORY ALLOCATION ERROR\n", me);
      break;
    default:
      printf("%d: INVALID RETURN CODE\n", me);
      break;
    }
  return;
}
int exGlobalSuccess(int rc, int nprocs, int me, int verbose)
{
  int fail = 0;
  int i;
  int *vals = (int *)malloc(nprocs * sizeof(int));

  MPI_Allgather(&rc, 1, MPI_INT, vals, 1, MPI_INT, MPI_COMM_WORLD);

  for (i=0; i<nprocs; i++){
    if (vals[i] != ZOLTAN_OK){
      if (verbose){
        exShowError(vals[i], "Result on process ", i);
      }
      fail = 1;
    }
  }

  free(vals);
  return fail;
}
void exPrintGlobalResult(char *s, int nprocs, int me,
              int begin, int import, int export, int change)
{
  int i;
  int *v1 = (int *)malloc(4 * sizeof(int));
  int *v2 = NULL;
  int *v;

  v1[0] = begin;
  v1[1] = import;
  v1[2] = export;
  v1[3] = change;

  if (me == 0){
    v2 = (int *)malloc(4 * nprocs * sizeof(int));
  }

  MPI_Gather(v1, 4, MPI_INT, v2, 4, MPI_INT, 0, MPI_COMM_WORLD);

  if (me == 0){
    fprintf(stdout,"======%s======\n",s);
    for (i=0, v=v2; i<nprocs; i++, v+=4){
      fprintf(stdout,"%d: originally had %d, import %d, export %d, %s\n",
        i, v[0], v[1], v[2],
        v[3] ? "a change of partitioning" : "no change");
    }
    fprintf(stdout,"==========================================\n");
    fflush(stdout); 

    free(v2);
  }

  free(v1);
}

#ifdef __cplusplus
} /* closing bracket for extern "C" */
#endif
