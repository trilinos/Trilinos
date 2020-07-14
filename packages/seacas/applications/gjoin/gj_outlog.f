C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C -*- Mode: fortran -*-
C=======================================================================
C $Id: outlog.f,v 1.1 1999/01/18 19:21:24 gdsjaar Exp $
C $Log: outlog.f,v $
C Revision 1.1  1999/01/18 19:21:24  gdsjaar
C ExodusII version of gjoin, needs testing and syncing with exodus 1 version, but is being committed to permit easier testing and modifications.  This was created by Dave Fry at Goodyear
C
c Revision 1.1.1.1  1998/11/05  16:23:27  a294617
c Initial import == gjoin 1.36
c
C Revision 1.1.1.1  1990/11/12 14:35:27  gdsjaar
C GJOIN - X1.00.40 - 7/17/90
C
c Revision 1.1  90/11/12  14:35:26  gdsjaar
c Initial revision
c
      SUBROUTINE OUTLOG (KLOG, NUMFLD, INTYP, CFIELD, IFIELD, RFIELD)
C=======================================================================

      CHARACTER*(*) CFIELD(*)
      INTEGER       IFIELD(*), INTYP(*)
      REAL          RFIELD(*)
      CHARACTER*132 STRING

      IF (KLOG .LE. 0) RETURN
      STRING = ' '

      DO 10 IFLD = 1, NUMFLD
         IF (INTYP(IFLD) .LT. 0) THEN
            CALL FFADDC (' ', STRING)
         ELSE IF (INTYP(IFLD) .EQ. 0) THEN
            CALL FFADDC (CFIELD(IFLD), STRING)
         ELSE IF (INTYP(IFLD) .EQ. 1) THEN
            CALL FFADDR (RFIELD(IFLD), STRING)
         ELSE IF (INTYP(IFLD) .EQ. 2) THEN
            CALL FFADDI (IFIELD(IFLD), STRING)
         ELSE
            CALL PRTERR ('PROGRAM', 'Unrecognized field type in OUTLOG')
         END IF
   10 CONTINUE

      WRITE (KLOG, '(A)') STRING(:LENSTR(STRING))

      RETURN
      END
