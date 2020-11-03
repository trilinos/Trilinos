C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBISTE (NDB, OPTION, ISTEP,
     &  NVARGL,
     *  NVARNP, NUMNP,
     *  NVAREL, NELBLK, NUMELB, ISEVOK, IDELB,
     *  NVARNS, NUMNPS, NNNPS,  ISNSOK, IDNPS,
     *  NVARSS, NUMESS, NEESS,  ISSSOK, IDESS,
     *  TIME, VARGL, VARNP, VAREL, VARNS, VARSS, *)
C=======================================================================
C     --*** DBISTE *** (EXOLIB) Read database variables for one time step
C     --   Written by Amy Gilkey - revised 10/14/87
C     --
C     --DBISTE reads the database global, nodal, and element variables
C     --for one time step.
C     --
C     --Parameters:
C     --   NDB - IN - the database number
C     --   OPTION - IN - ' ' to not store, '*' to store all, else store options:
C     --      'G' to store global variables
C     --      'E' to store element variables
C     --      'N' to store nodal variables
C     --      'M' to store nodeset variables
C     --      'S' to store sideset variables
C     --   ISTEP - IN - the time step number
C     --   NVARGL - IN - the number of global variables
C     --   NVARNP - IN - the number of nodal variables
C     --   NUMNP - IN - the number of nodes
C     --   NVAREL - IN - the number of element variables
C     --   NELBLK - IN - the number of element blocks
C     --   NUMELB - IN - the number of elements per block
C     --   ISEVOK - IN - the element block variable truth table;
C     --      variable i of block j exists iff ISEVOK(j,i)
C     --   TIME - OUT - the time step time
C     --   VARGL - OUT - the global variables for the time step (if OPTION)
C     --   VARNP - OUT - the nodal variables for the time step (if OPTION)
C     --   VAREL - OUT - the element variables for the time step (if OPTION)
C     --   * - return statement if error encountered, including end-of-file;
C     --      message is printed
C     --
C     --Database must be positioned in front of time step upon entry;
C     --upon exit positioned after time step.

      include 'exodusII.inc'
      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER ISTEP
      INTEGER NVARGL
      INTEGER NVARNP, NUMNP
      INTEGER NVAREL, NELBLK, NUMELB(*), IDELB(*), ISEVOK(*)
      INTEGER NVARNS, NUMNPS, NNNPS(*),  IDNPS(*), ISNSOK(*)
      INTEGER NVARSS, NUMESS, NEESS(*),  IDESS(*), ISSSOK(*)
      REAL TIME
      REAL VARGL(*)
      REAL VARNP(*)
      REAL VAREL(*)
      REAL VARNS(*)
      REAL VARSS(*)

      CHARACTER*80 ERRMSG

C     --Read step time
      CALL EXGTIM (NDB, ISTEP, TIME, IERR)

C     --Read global variables

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'G') .GT. 0)) THEN
         if (nvargl .gt. 0) then
           call exggv (ndb, istep, nvargl, vargl, ierr)
           if (ierr .lt. 0) goto 180
         end if
      END IF

C     --Read nodal variables

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'N') .GT. 0)) THEN
         ioff = 1
         do 150 ivar = 1, nvarnp
           call exgnv (ndb, istep, ivar, numnp, varnp(ioff), ierr)
           if (ierr .lt. 0) goto 190
           ioff = ioff + numnp
 150     continue
      END IF

C     --Read element variables
      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'E') .GT. 0)) THEN
         CALL DBIST2 (NDB, 'E', ISTEP, NVAREL, MAX(NVAREL,1), NELBLK,
     &        MAX(NELBLK,1), ISEVOK, VAREL, NUMELB,
     &        IDELB, IVAR, IELB, *200)
      END IF

C     --Read nodeset variables
      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'M') .GT. 0)) THEN
         CALL DBIST2 (NDB, 'M', ISTEP, NVARNS, MAX(NVARNS,1), NUMNPS,
     &        MAX(NUMNPS,1), ISNSOK, VARNS, NNNPS,
     &        IDNPS, IVAR, IELB, *210)
      END IF

C     --Read sideset variables
      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'S') .GT. 0)) THEN
         CALL DBIST2 (NDB, 'S', ISTEP, NVARSS, MAX(NVARSS,1), NUMESS,
     &        MAX(NUMESS,1), ISSSOK, VARSS, NEESS,
     &        IDESS, IVAR, IELB, *220)
      END IF

      RETURN

 180  CONTINUE
      WRITE (ERRMSG, '(A, I8)', IOSTAT=IDUM)
     &     'GLOBAL VARIABLES for TIME STEP', ISTEP
      GOTO 240

 190  CONTINUE
      WRITE (ERRMSG, '(A, I8, A, I8)', IOSTAT=IDUM)
     &     'NODAL VARIABLE', IVAR, ' for TIME STEP', ISTEP
      GOTO 240

 200  CONTINUE
      WRITE (ERRMSG, '(A, I8, A, I8, A, I8)', IOSTAT=IDUM)
     &     'ELEMENT VARIABLE', IVAR, ' of BLOCK', IELB,
     &     ' for TIME STEP', ISTEP
      GOTO 240

 210  CONTINUE
      WRITE (ERRMSG, '(A, I8, A, I8, A, I8)', IOSTAT=IDUM)
     &     'NODESET VARIABLE', IVAR, ' of NODESET', IELB,
     &     ' for TIME STEP', ISTEP
      GOTO 240

 220  CONTINUE
      WRITE (ERRMSG, '(A, I8, A, I8, A, I8)', IOSTAT=IDUM)
     &     'SIDESET VARIABLE', IVAR, ' of SIDESET', IELB,
     &     ' for TIME STEP', ISTEP
      GOTO 240

 240  CONTINUE
      CALL DBERR (IERR, ERRMSG)
      RETURN 1
      END
