/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/* cgidmy - dummy routine loaded if no load line devices specified
 * 28 Nov 1988 - last date modified
 * @(#) cgidmy.c version 1.1 of 2/24/89 10:20:59
 * Pat McGee, jpm@lanl.gov
 */
/* 4/10/89 - changed name to cgisdum.c */

#include <stdio.h>
#include <stdlib.h>

void cgi_def_ini(void)
{
  fprintf(stderr, "cgi: error: no device drivers specified\n");
  exit(1);
} /* end cgi_def_ini */
