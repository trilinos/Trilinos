/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

#include "smalloc.h" // for smalloc
#include "structs.h" // for orthlink, orthlink_float

/* Allocate space for new orthlink, double version. */
struct orthlink *makeorthlnk(void)
{
  struct orthlink *newlnk;

  newlnk = smalloc(sizeof(struct orthlink));
  return (newlnk);
}

/* Allocate space for new orthlink, float version. */
struct orthlink_float *makeorthlnk_float(void)
{
  struct orthlink_float *newlnk;

  newlnk = smalloc(sizeof(struct orthlink_float));
  return (newlnk);
}
