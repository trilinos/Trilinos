#ifndef CHACO_INTERNAL_INTERNAL_H
#define CHACO_INTERNAL_INTERNAL_H

/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/* Bidirectional list of integers. */
struct bidint
{
  int            val;
  struct bidint *prev;
  struct bidint *next;
};

#endif
