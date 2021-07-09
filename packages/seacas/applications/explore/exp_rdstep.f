C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE RDSTEP (NDB, NCSTEP,
     &     NUMELB, IDELB, ISEVOK,
     &     TIME, VARGL, VARNP, VAREL,
     $     neldm)
C=======================================================================

C   --*** RDSTEP *** (EXPLORE) Read current database variables
C   --
C   --RDSTEP reads the variable data for the next time step from the
C   --database.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   NCSTEP - IN/OUT - the current step number
C   --   NVARGL - IN - the number of global variables
C   --   NVARNP - IN - the number of nodal variables
C   --   NVAREL - IN - the number of element variables
C   --   NUMELB - IN - the number of elements per block
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(i,j) is NOT 0
C   --   TIME - OUT - the time step time
C   --   VARGL - OUT - the global variables for the time step
C   --   VARNP - OUT - the nodal variables for the time step
C   --   VAREL - OUT - the element variables for the time step

      include 'exodusII.inc'
      include 'exp_dbnums.blk'

      INTEGER NUMELB(*)
      INTEGER IDELB(*)
      INTEGER ISEVOK(nvarel,*)
      REAL VARGL(*)
      REAL VARNP(numnp, *)
      REAL VAREL(numel, *)

      CHARACTER*80 ERRMSG

C ... Initialize variables to known value in case error on read
      TIME = -999.9
      CALL INIREA (NVARGL, -999.9, VARGL)
      CALL INIREA (NVARNP * NUMNP, -999.9, VARNP)
      CALL INIREA (NVAREL * NUMEL, -999.9, VAREL)

C ... Read the time for the specified step
      call exgtim(ndb, ncstep, time, ierr)
      if (ierr .ne. 0) go to 150

C ... Read the global variables (if any)
      if (nvargl .gt. 0) then
         call exggv (ndb, ncstep, nvargl, vargl, ierr)
         if (ierr .ne. 0) go to 170
      end if

C ... Read the nodal variables (if any)
      DO 120 I = 1, NVARNP
        call exgnv (ndb, ncstep, i, numnp, varnp(1,i), ierr)
        if (ierr .ne. 0) go to 180
 120  CONTINUE

C ... Read the element variables (if any)
      DO 140 I = 1, NVAREL
        IEL0 = 1
        DO 130 IELB = 1, NELBLK
          IF (ISEVOK(I,IELB) .NE. 0 .AND. numelb(ielb) .GT. 0) THEN
             call exgev (ndb, ncstep, i, idelb(ielb),
     &            numelb(ielb), varel(iel0,i), ierr)
             if (ierr .ne. 0) go to 190
          END IF
          IEL0 = IEL0 + NUMELB(IELB)
 130    CONTINUE
 140  CONTINUE

      RETURN

  150 CONTINUE
      WRITE (ERRMSG, 10000, IOSTAT=IDUM)
     &   'TIME on step', NCSTEP
      GOTO 200
  170 CONTINUE
      WRITE (ERRMSG, 10000, IOSTAT=IDUM)
     &   'GLOBAL VARIABLES on step', NCSTEP
      GOTO 200
  180 CONTINUE
      WRITE (ERRMSG, 10000, IOSTAT=IDUM)
     &   'NODAL VARIABLE', I, ' on step', NCSTEP
      GOTO 200
  190 CONTINUE
      WRITE (ERRMSG, 10000, IOSTAT=IDUM)
     &   'ELEMENT VARIABLE', I, ' in block', IELB, ' on step', NCSTEP
      GOTO 200
  200 CONTINUE
      CALL WDBERR (IERR, ERRMSG)

      RETURN

10000  FORMAT (5 (A, I12))
      END
