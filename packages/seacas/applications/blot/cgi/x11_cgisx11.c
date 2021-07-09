/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/* cgisx11.c - linker specifiable driver routine for driver
 *   X.11 (x11)
 * Sandia National Laboratories, Div 2634
 * Sun Nov 19 12:02:52 MST 1989 - last date modified
 */

#include "ifdefx.h"
#include "mdcgi.h"

void cgix11_(); /* tell linker to load driver */

void cgisx11() /* make name external so linker will load file*/ {}

void cgi_def_ini()
{
  anything *devid;

  xcact_(cgix11_, &devid);
  xcoon_(&devid);

} /* end cgi_def_ini */
