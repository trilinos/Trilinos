/*****************************************************************************
 * Zoltan Library for Parallel Applications                                  *
 * Copyright (c) 2008, Sandia National Laboratories.                         *
 * For more info, see the README file in the top-level Zoltan directory.     *
 *****************************************************************************/
/*****************************************************************************
 * CVS File Information :
 *    $RCSfile$
 *    $Author$
 *    $Date$
 *    $Revision$
 ****************************************************************************/


#ifdef __cplusplus
/* if C++, define the rest of this header file as extern C */
extern "C" {
#endif

#define MM_INPUT_C

#include "common.h"
#include "scotch.h"
#include "mm_input.h"

/*****************************************************/
/*                                                   */
/* These routines handle source matrix market Input. */
/*                                                   */
/*****************************************************/


int
readMMHeader(
	     SCOTCH_Num * vertnbr,
	     SCOTCH_Num * edgenbr,
FILE * const stream)
{
  int c;
  SCOTCH_Num vertnum0, vertnum1;
  int symmetric = 0;
  char firstline[1024];
  int pos = 0;

  memSet(firstline, '\0', 1024);
  fgets(firstline, 16, stream);

  if (strcmp(firstline, "%%MatrixMarket "))
    return (-1);

  while (((c = fgetc(stream)) != '\n') && (c != EOF)) {
    firstline[pos++] = toupper(c);
    if (pos == 1024)
      return (-1);
  }
  firstline[pos] = '\0';

  if (strstr(firstline, "SYMMETRIC"))
    symmetric = 1;

  while ((c = fgetc(stream)) == '%') {
    if (skipLine(stream) != 0)
      return (-1);                                /* End of file reached */
  }
  ungetc (c, stream);

  if (intLoad (stream, &vertnum0) != 1) {         /* Read row number */
    return (-1);
  }
  if (intLoad (stream, &vertnum1) != 1) {         /* Read col number */
    return (-1);
  }
  if (vertnum0 != vertnum1) {                     /* Non square matrix */
      return (-1);
  }
  *vertnbr = vertnum1;
  if (intLoad (stream, edgenbr) != 1) { /* Read edge number */
    return (-1);
  }
  *edgenbr -= *vertnbr;           /* No loop in graph */
  return (symmetric);
}

int skipLine(FILE * const stream)
{
  int c;

  while (((c = fgetc(stream)) != '\n') && (c != EOF));

  if (c == EOF)
    return (1);
  return (0);
}


/* This routine loads a source graph from
** the given stream, corresponding to a MatrixMarket file.
** It returns:
** - 0   : on success.
** - !0  : on error.
*/

/**
 *  Note:
 *  The algorithm can be more efficient in term of memory and time
 *  for the symmetric graphs
 */

int
graphLoadMM (
SCOTCH_Graph * restrict const      grafptr,              /* Graph structure to fill              */
FILE * const                stream)               /* Stream from which to read graph data */
{
  SCOTCH_Num          baseval;
  SCOTCH_Num *        verttax;
  SCOTCH_Num *        arctax;
  GraphEdge *         edgetax;
  SCOTCH_Num          edgenum;                    /* Number of edges really allocated */
  SCOTCH_Num          vertnum;
  SCOTCH_Num          vertnbr;
  SCOTCH_Num          edgenbr;
  int                 symmetry = 0;

  SCOTCH_Num          degree;
  SCOTCH_Num          neighbor;
  SCOTCH_Num          edgeidx;

  baseval = 1;                                    /* Matrix starts at 1 1 */

  if ((symmetry = readMMHeader(&vertnbr, &edgenbr, stream)) < 0) {
    errorPrint ("graphLoadMM: Invalid input file");
    return (1);
  }

  if ((verttax  = (SCOTCH_Num *)memAlloc((size_t) ((vertnbr + 1) * sizeof (SCOTCH_Num)))) == NULL){
    errorPrint ("graphLoadMM: out of memory (1)");
    return     (1);
  }
  verttax -= baseval;

  if ((edgetax = (GraphEdge*) memAlloc ((size_t) (2 * edgenbr * sizeof(GraphEdge)))) == NULL) {
    memFree (verttax);
    errorPrint ("graphLoadMM: out of memory (2)");
    return     (1);
  }

  edgetax -= baseval ;


  for (edgenum = baseval; ; ) {          /* Read the file and fill an array of edges */
    SCOTCH_Num vertnum0, vertnum1;

    if (intLoad (stream, &vertnum0) != 1) {        /* Read vertex degree */
      break;
    }
    if (intLoad (stream, &vertnum1) != 1) {        /* Read vertex degree */
      break;
    }

    skipLine (stream);
    if (vertnum0 == vertnum1)
      continue;

    edgetax[edgenum].vert[0] = vertnum0;
    edgetax[edgenum].vert[1] = vertnum1;
    edgenum++;

    edgetax[edgenum].vert[0] = vertnum1;          /* Perhaps only (vertnum0, vertnum1) is defined */
    edgetax[edgenum].vert[1] = vertnum0;
    edgenum++;
  }
  edgenbr = (edgenum - baseval);

  intSort2asc2 (edgetax + baseval, edgenbr);

  if ((edgetax = (GraphEdge *)memRealloc (edgetax + baseval,
					  (size_t) (edgenbr * sizeof(GraphEdge)))) == NULL) {
    memFree (verttax);
    errorPrint ("graphLoadMM: out of memory (3)");
    return     (1);
  }
  edgetax -= baseval;
  if ((arctax = (SCOTCH_Num*) memAlloc ((size_t) (edgenbr * sizeof(SCOTCH_Num)))) == NULL) {
    memFree (edgetax + baseval);
    memFree (verttax + baseval);
    errorPrint ("graphLoadMM: out of memory (4)");
    return     (1);
  }
  arctax -= baseval;

  for (edgeidx = edgenum = baseval, vertnum=baseval - 1, neighbor=-1, degree = 0 ;
       edgeidx < edgenbr + baseval ; ++edgeidx) {
    if (edgetax[edgeidx].vert[0] != vertnum) {    /* vertex not yet processed */
      for (++vertnum; vertnum <  edgetax[edgeidx].vert[0] ; ++vertnum)
	verttax[vertnum] = edgenum;      /* Take care of only diagonal term in column */
      verttax[vertnum] = edgenum;
      neighbor = -1;
    }
    if (neighbor != edgetax[edgeidx].vert[1]) {   /* Edge not yet processed */
      degree++;
      neighbor = arctax[edgenum++] = edgetax[edgeidx].vert[1];
    }
  }
  memFree(edgetax + baseval);
  edgenbr = edgenum - baseval;

  if ((arctax = (SCOTCH_Num*) memRealloc (arctax + baseval,
					  (size_t) (edgenbr * sizeof(SCOTCH_Num)))) == NULL) {
    memFree (verttax + baseval);
    errorPrint ("graphLoadMM: out of memory (5)");
    return     (1);
  }
  arctax -= baseval;

  return (SCOTCH_graphBuild (grafptr, baseval, vertnbr, verttax + baseval,
			     NULL, NULL, NULL, edgenbr, arctax + baseval, NULL));
}

#ifdef __cplusplus
}
#endif
