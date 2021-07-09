/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/* cgispst.c - linker specifiable driver routine for driver Postscript B/W (pst) */
#include "ifdefx.h"
#include "mdcgi.h"  // for xcact_, xcoon_
#include "stdtyp.h" // for anything

void cgipst_(); /* tell linker to load driver */

void cgispst(void) /* make name external so linker will load file*/ {}

void cgi_def_ini(void)
{
  anything *devid;

  xcact_(cgipst_, &devid);
  xcoon_(&devid);

} /* end cgi_def_ini */
