C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PRHIST (OPTION, NOUT, NVARHI, LISHV, NAMEHV, VARHI)
C=======================================================================

C   --*** PRHIST *** (BLOT) Display current database history variables
C   --   Written by Amy Gilkey - revised 11/05/87
C   --
C   --PRHIST displays the history data for a time step.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options
C   --   NOUT - IN - the output file, <=0 for standard
C   --   NVARHI - IN - the number of history variables
C   --   LISHV - IN - the indices of the selected history variables
C   --   NAMEHV - IN - the names of the history variables
C   --   VARHI - IN - the history variables for the time step

      CHARACTER*(*) OPTION
      INTEGER LISHV(0:*)
      CHARACTER*(*) NAMEHV(*)
      REAL VARHI(*)

      IF (NOUT .GT. 0) WRITE (NOUT, 10000)

      IF (NOUT .GT. 0) THEN
         WRITE (NOUT, 10010) (NAMEHV(LISHV(I)), I=1,LISHV(0))
      ELSE
         WRITE (*, 10010) (NAMEHV(LISHV(I)), I=1,LISHV(0))
      END IF

      IF (NOUT .GT. 0) THEN
         WRITE (NOUT, 10020, IOSTAT=IDUM)
     &      (VARHI(LISHV(I)), I=1,LISHV(0))
      ELSE
         WRITE (*, 10020, IOSTAT=IDUM)
     &      (VARHI(LISHV(I)), I=1,LISHV(0))
      END IF

      RETURN

10000  FORMAT (/, 1X, 'HISTORY TIME STEP VARIABLES')
10010  FORMAT (/, 1X, 9X, 5 (3X, A, 3X), :, /,
     &   (1X, 9X, 5 (3X, A, 3X)))
10020  FORMAT (1X, 'History', 2X, 5 (1X, E13.6), :, /,
     &   (1X, 9X, 5 (1X, E13.6)))
      END
