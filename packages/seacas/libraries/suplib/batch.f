C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION BATCH ()
C=======================================================================
C$Id: batch.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: batch.f,v $
CRevision 1.2  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:12:03  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:12:02  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:03  gdsjaar
c Initial revision
c

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
