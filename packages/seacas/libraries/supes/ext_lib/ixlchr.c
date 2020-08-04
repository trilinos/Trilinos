/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */
/*
 */

#include "fortranc.h"
#if defined(ADDC_)
FTNINT ixlchr_(long *chrvar)
#else
FTNINT ixlchr(long *chrvar)
#endif
{
  return ((FTNINT)chrvar);
}

/*

************************************************************************
C     DESCRIPTION:
C     This function returns the absolute location of a character variable.
C     This location must be measured in character storage units.

C     FORMAL PARAMETERS:
C     CHRVAR    CHARACTER       Character Variable
************************************************************************

*/
