C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION BATCH ()
C=======================================================================

C   --*** BATCH *** (ETCLIB) Return batch versus interactive flag
C   --   Written by Amy Gilkey - revised 01/20/87
C   --
C   --BATCH returns true iff the program is in batch rather that interactive
C   --mode.

C   --Routines Called:
C   --   EXPARM - (SUPES) Get batch vs. interactive flag

      CHARACTER*8 CDUM

      LOGICAL FIRST, SVBATC
      SAVE FIRST, SVBATC

      DATA FIRST / .TRUE. /

      IF (FIRST) THEN
         CALL EXPARM (CDUM, CDUM, IMODE, IDUM, IDUM, IDUM)
         SVBATC = (IMODE .EQ. 0)
         FIRST = .FALSE.
      END IF

      BATCH = SVBATC

      RETURN
      END
