/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "smalloc.h"
#include <stdio.h>

static void assign_out_normal(int   nvtxs,  /* number of vertices to output */
                              int  *sets,   /* values to be printed */
                              char *outname /* name of output file */
)
{
  FILE *fout; /* output file */
  int   i;    /* loop counter */

  /* Print assignment in simple format. */

  if (outname != NULL) {
    fout = fopen(outname, "w");
  }
  else {
    fout = stdout;
  }

  for (i = 1; i <= nvtxs; i++) {
    fprintf(fout, "%d\n", sets[i]);
  }

  if (outname != NULL) {
    fclose(fout);
  }
}

static void assign_out_inv(int   nvtxs,  /* number of vertices to output */
                           int  *sets,   /* values to be printed */
                           int   nsets,  /* number of sets */
                           char *outname /* name of output file */
)
{
  FILE *fout;    /* output file */
  int  *size;    /* # vtxs in sets / index into inorder */
  int  *inorder; /* list of vtxs in each set */
  int   i, j;    /* loop counter */

  /* Print assignment in inverted format. */
  /* For each set, print # vertices, followed by list. */

  if (outname != NULL) {
    fout = fopen(outname, "w");
  }
  else {
    fout = stdout;
  }

  size    = smalloc((nsets + 1) * sizeof(int));
  inorder = smalloc(nvtxs * sizeof(int));

  for (j = 0; j < nsets; j++) {
    size[j] = 0;
  }
  for (i = 1; i <= nvtxs; i++) {
    ++size[sets[i]];
  }

  /* Modify size to become index into vertex list. */
  for (j = 1; j < nsets; j++) {
    size[j] += size[j - 1];
  }
  for (j = nsets - 1; j > 0; j--) {
    size[j] = size[j - 1];
  }
  size[0] = 0;

  for (i = 1; i <= nvtxs; i++) {
    j                = sets[i];
    inorder[size[j]] = i;
    ++size[j];
  }

  /* The inorder array now clumps all the vertices in each set. */
  /* Now reconstruct size array to index into inorder. */
  for (j = nsets - 1; j > 0; j--) {
    size[j] = size[j - 1];
  }
  size[0]     = 0;
  size[nsets] = nvtxs;

  /* Now print out the sets. */
  for (j = 0; j < nsets; j++) {
    fprintf(fout, "  %d\n", size[j + 1] - size[j]);
    for (i = size[j]; i < size[j + 1]; i++) {
      fprintf(fout, "%d\n", inorder[i]);
    }
  }

  if (outname != NULL) {
    fclose(fout);
  }
  sfree(size);
  sfree(inorder);
}

void assign_out(int   nvtxs,  /* number of vertices to output */
                int  *sets,   /* values to be printed */
                int   nsets,  /* number of sets */
                char *outname /* name of output file */
)
{
  extern int OUT_ASSIGN_INV; /* print assignment in inverted form? */

  if (OUT_ASSIGN_INV) {
    assign_out_inv(nvtxs, sets, nsets, outname);
  }
  else {
    assign_out_normal(nvtxs, sets, outname);
  }
}
