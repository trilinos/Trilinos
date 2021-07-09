/*
 * Copyright(C) 1999-2020 National Technology & Engineering Solutions
 * of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
 * NTESS, the U.S. Government retains certain rights in this software.
 *
 * See packages/seacas/LICENSE for details
 */

/*
 *   This routine is here for the purpose of assigning a replacement
 *   file name thus overriding the default, machine dependent, FORTRAN
 *   I/O file assignments.  For example, on the VAX, under VAX
 *   FORTRAN, an OPEN(n) statement writes the data to the file
 *   FOR00n.DAT in the default directory.  One can override this
 *   filename by ASSIGN'ing the logical name FOR0NN to a filename
 *   before executing the image.  Sadly, the CRAY version of FORTRAN
 *   does not allow for such an eventuality.  So, to affect the
 *   changes, one must manually check for the file name via this
 *   routine which tests for the existence of a execution line string
 *   of the form "tapen=filename".  Furthermore, one must test the
 *   return filename length to see if such a string does exist and
 *   loop to the appropriate form of OPEN() statement.  (Note that for
 *   the time being, this capability under CTSS is *NOT* part of this
 *   module.  Maybe it can be fixed later.)  For example, consider the
 *   following code segment:
 *
 *        CALL EXNAME( 7, FILENM, LN )
 *        IF( LN .EQ. 0 ) THEN             ! Then just take the default
 *            OPEN( 7 )
 *        ELSE                             ! I've found a file name, use it.
 *            OPEN( 7, FILE=FILENM )
 *        ENDIF
 *
 *   In fact, the same sort of thing happens in the UNIX world.
 *   To provide as least some compatibility with the VAX world,  we test
 *   for the ENVIRONMENT variables of the form FOR0NN.
 *
 */
#include "fortranc.h"
#include <stdio.h>
#include <stdlib.h> /* getenv */
#include <string.h> /* strlen */

#define TRUE 1
#define FALSE 0

#define INSYMBOL "FOR0"
#define EXSYMBOL "EXT"

static char *copy_string(char *dest, char const *source, long int elements)
{
  char *d;
  for (d = dest; d + 1 < dest + elements && *source; d++, source++) {
    *d = *source;
  }
  *d = '\0';
  return d;
}

#if defined(ADDC_)
void exname_(FTNINT *iunit, char *name, FTNINT *ln, long int nlen)
#else
void exname(FTNINT *iunit, char *name, FTNINT *ln, long int nlen)
#endif
{
  char string[3];
  char symbol[7];
  int  ExtSymbol;

  char *darg;

  if (*iunit > -100 && *iunit < 100) {
#if Build64
    sprintf(string, "%02ld", labs(*iunit));
#else
    sprintf(string, "%02d", abs(*iunit));
#endif

    if (*iunit > 0) {
      ExtSymbol = FALSE;
      sprintf(symbol, "%s%s", INSYMBOL, string);
    }
    else {
      ExtSymbol = TRUE;
      sprintf(symbol, "%s%s", EXSYMBOL, string);
    }

    if ((darg = (char *)getenv(symbol)) != (char *)NULL) {
      unsigned int DargLen = strlen(darg);
      /* We need this to give us the length of the ENVIRONMENT
       * variable while calling strlen() only once.
       */
      copy_string(name, darg, DargLen + 1);
      *ln = DargLen;
    }
    else if (!ExtSymbol) {
#if Build64
      sprintf(name, "fort.%ld", labs(*iunit));
#else
      sprintf(name, "fort.%d", abs(*iunit));
#endif
      *ln = strlen(name);
    }
    else {
      /* Then I have referenced an external symbol that has not been
       *  defined...
       */
      *name = '\0';
      *ln   = 0;
    }
  }
  else {
#if Build64
    fprintf(stderr,
            "ERROR: SUPES exname - invalid unit number %ld.  Must be between -100 and 100.\n",
            *iunit);
#else
    fprintf(stderr,
            "ERROR: SUPES exname - invalid unit number %d.  Must be between -100 and 100.\n",
            *iunit);
#endif
  }
}
/*
************************************************************************
C     DESCRIPTION:
C     This routine returns either a valid file name for a logical unit
C     number or the value of a symbol. If IUNIT .GT. 0, the returned
C     name may be used as the value for the FILE specification in an
C     OPEN statement. If IUNIT .LE. 0, the returned name is the value of
C     a system symbol. It is assumed that the unit/file name or symbol
C     linkage will be passed to this routine during program execution.
C     A null string (LN = 0) will be returned if no name is available.

C     FORMAL PARAMETERS:
C     IUNIT     INTEGER         Logical Unit Number ( >0 )
C                            or Symbol ID ( = -IUNIT )
C     NAME      CHARACTER       File Name
C     LN        INTEGER         Length of File Name

************************************************************************
*/
