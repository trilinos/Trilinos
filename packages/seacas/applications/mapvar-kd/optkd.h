/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
#ifndef OPTKD_H
#define OPTKD_H
#define BUCKETSIZE 100

typedef struct optkdnode
{
  int               bucket;  /* true if the node is a bucket node */
  int               discrim; /* discriminator of node */
  real              cutval;
  struct optkdnode *loson, *hison;
  int               lopt, hipt;
} optkdNode;
#endif
