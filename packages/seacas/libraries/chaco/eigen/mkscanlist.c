/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "smalloc.h" // for smalloc
#include "structs.h" // for scanlink
#include <stdio.h>   // for NULL

struct scanlink *mkscanlist(int depth)
{
  struct scanlink *prevlnk;
  struct scanlink *newlnk;
  int              i;

  prevlnk       = smalloc(sizeof(struct scanlink));
  prevlnk->pntr = NULL;
  newlnk        = prevlnk; /* in case the list is one long */
  for (i = 1; i <= (depth - 1); i++) {
    newlnk       = smalloc(sizeof(struct scanlink));
    newlnk->pntr = prevlnk;
    prevlnk      = newlnk;
  }
  return (newlnk);
}
