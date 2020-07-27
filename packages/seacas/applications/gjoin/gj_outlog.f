C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C -*- Mode: fortran -*-
C=======================================================================
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
