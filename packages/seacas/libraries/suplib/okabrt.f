C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION OKABRT (ISOK)
C=======================================================================

C   --*** OKABRT *** (ETCLIB) Initialize cancel function
C   --   Written by Amy Gilkey - revised 12/21/87
C   --
C   --OKABRT initializes the cancel flag.  It must be called before ISABRT.

C   --Routines Called:
C   --   CPUIFC - (PLTLIB) Check interrupt flag

      LOGICAL ISABRT
      LOGICAL ISOK

      LOGICAL CPUIFC, LDUM

      LOGICAL DOABRT
      SAVE DOABRT

      DATA DOABRT / .FALSE. /

C   --Initialize enable cancel flag
      DOABRT = ISOK

      IF (DOABRT) THEN
C      --Initialize cancel flag
         LDUM = CPUIFC (.TRUE.)
      END IF

      OKABRT = DOABRT

      RETURN

C=======================================================================
      ENTRY ISABRT ()
C=======================================================================

C   --*** ISABRT *** (ETCLIB) Check cancel function
C   --   Written by Amy Gilkey - revised 12/17/87
C   --
C   --ISABRT checks the cancel flag.  If it is set, it aborts the current
C   --processing.  In any case, the value of the cancel flag is returned
C   --as the function value.

C   --Routines Called:
C   --   CPUIFC - (PLTLIB) Check interrupt flag

      IF (DOABRT) THEN
C      --Return cancel flag
         ISABRT = CPUIFC (.FALSE.)

         IF (ISABRT) THEN
C         --Send abort message
            WRITE (*, '(1X, A)') '*** Processing aborted ***'
         END IF

      ELSE
         ISABRT = .FALSE.
      END IF

      RETURN
      END
