/*
 * Copyright(C) 1999-2020, 2022 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/* Bidirectional list of integers. */
#pragma once

struct bidint
{
  int            val;
  struct bidint *prev;
  struct bidint *next;
};
