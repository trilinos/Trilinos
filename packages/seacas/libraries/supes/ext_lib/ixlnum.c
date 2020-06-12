/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 * 
 * See packages/seacas/LICENSE for details
 */
/*
 * $Id: ixlnum.c,v 1.14 2008/03/14 13:22:39 gdsjaar Exp $
 */

#include "fortranc.h"

#if defined(ADDC_)
FTNINT ixlnum_(long *numvar)
#else
FTNINT ixlnum(long *numvar)
#endif
{
  return ((long)numvar / sizeof(FTNREAL));
}

/*
      INTEGER FUNCTION IXLNUM( NUMVAR )
C
************************************************************************
C
C     DESCRIPTION:
C     This function returns the absolute location of a numeric variable.
C     This location must be measured in numeric storage units.
C
C     FORMAL PARAMETERS:
C     NUMVAR    INTEGER         Numeric Variable
C
************************************************************************
*/
