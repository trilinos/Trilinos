C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE TOSTEP (NSTEP, NUMELB, IDELB, ISEVOK,
     &   TIME, VARGL, VARNP, VAREL)
C=======================================================================

C   --*** TOSTEP *** (EXPLORE) Move database to given time step number
C   --
C   --TOSTEP reads the given time step variables (NSTEP), setting
C   --the time step to the new current time step.  A message is
C   --displayed if the time step is too large.
C   --
C   --Parameters:
C   --   NSTEP - IN/OUT - the number of the time step to be read;
C   --      the returned number of the time step actually read
C   --   NUMELB - IN - the number of elements per block
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(i,j) is NOT 0
C   --   TIME - OUT - the time step time
C   --   VARGL - IN/OUT - the global variables for the time step
C   --   VARNP - IN/OUT - the nodal variables for the time step
C   --   VAREL - IN/OUT - the element variables for the time step
C   --
C   --Common Variables:
C   --   Sets and uses NCSTEP of /DBASE/
C   --   Uses NELBLK, NQAREC, NINFO, NVARHI, NVARGL, NVARNP, NVAREL, NSTEPS
C   --      of /DBNUMS/

      include 'exp_dbase.blk'
      include 'exp_dbnums.blk'

      INTEGER NUMELB(*)
      INTEGER IDELB(*)
      INTEGER ISEVOK(NVAREL,*)
      REAL VARGL(*), VARNP(*), VAREL(*)

      CHARACTER*5 STRA

      LOGICAL FIRST
      SAVE FIRST

      INTEGER NCSTEP
      REAL CTIME
      save ctime
C      --NCSTEP - the current step number
C      --CTIME - the current time

      DATA FIRST / .TRUE. /

      IF (FIRST) THEN
         NCSTEP = -999
         FIRST = .FALSE.
      END IF

      IF (NSTEP .LE. 0) THEN
         CALL INTSTR (1, 0, NSTEP, STRA, LSTRA)
         CALL PRTERR ('CMDERR',
     &      'Invalid time step number ' // STRA(:LSTRA))
         GOTO 150
      END IF

      if (nstep .ne. ncstep) then
         CALL RDSTEP (NDB, NSTEP,
     &        NUMELB, IDELB, ISEVOK,
     &        CTIME, VARGL, VARNP, VAREL,
     $        max(1,nvarel))
         ncstep = nstep
      END IF
      TIME = CTIME

      IF (NSTEP .GT. NSTEPS) THEN
         IF (NCSTEP .EQ. 0) THEN
            WRITE (*, 10010) 'There are no time steps'
         ELSE
            WRITE (*, 10010)
     &         'All time steps have been read from the database'
         END IF
      END IF

  150 CONTINUE
      NSTEP = NCSTEP
      RETURN

10010  FORMAT (1X, 5A)
      END
