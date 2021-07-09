C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE TPCHK (PRTOK, OKAY)
C=======================================================================

C   --*** TPCHK *** (TPLOT) Check that TPLOT is a valid program to run
C   --   Written by Amy Gilkey - revised 11/03/87
C   --
C   --TPCHK checks that TPLOT is a valid program.  This is true if
C   --there are history, global, nodal or element variables in the database.
C   --
C   --If the program is valid, a welcome message is printed.
C   --
C   --Parameters:
C   --   PRTOK - IN - if true, print error messages and welcome messages,
C   --      if false, check but do not print anything
C   --   OKAY - OUT - true iff the program is valid
C   --
C   --Common Variables:
C   --   Uses NVARHI, NVARGL, NVARNP, NVAREL, NSTEPS, NSTEPW of /DBNUMS/

      include 'dbnums.blk'

      LOGICAL PRTOK
      LOGICAL OKAY

C   --Check header information from database

      IF ((NVARHI + NVARGL + NVARNP + NVAREL) .LE. 0) THEN
         IF (PRTOK) CALL PRTERR ('CMDERR',
     &      'No variables are defined')
         GOTO 100
      END IF
      IF ((NSTEPS .LE. 0) .OR.
     &   ((NSTEPW .LE. 0) .AND. (NVARHI .LE. 0))) THEN
         IF (PRTOK) CALL PRTERR ('CMDERR',
     &      'No time steps with variables are defined')
         GOTO 100
      END IF

C   --Write greeting message
      IF (PRTOK) WRITE (*, 10000)

      OKAY = .TRUE.
      RETURN

  100 CONTINUE
      OKAY = .FALSE.
      RETURN
10000  FORMAT (' TPLOT - a time history or X-Y plot program')
      END
