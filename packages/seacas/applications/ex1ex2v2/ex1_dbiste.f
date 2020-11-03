C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBISTE (NDB, OPTION, ISTEP,
     &   NVARHI, NVARGL, NVARNP, NUMNP, NVAREL, NELBLK, NUMELB, ISEVOK,
     &   TIME, WHOTIM, VARHI, VARGL, VARNP, VAREL, *)
C=======================================================================
C   --*** DBISTE *** (EXOLIB) Read database variables for one time step
C   --
C   --DBISTE reads the database global, nodal, and element variables
C   --for one time step.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   OPTION - IN - ' ' to not store, '*' to store all, else store options:
C   --      'H' to store history variables
C   --      'G' to store global variables
C   --      'E' to store element variables
C   --      'N' to store nodal variables
C   --   ISTEP - IN - the time step number
C   --   NVARHI - IN - the number of history variables
C   --   NVARGL - IN - the number of global variables
C   --   NVARNP - IN - the number of nodal variables
C   --   NUMNP - IN - the number of nodes
C   --   NVAREL - IN - the number of element variables
C   --   NELBLK - IN - the number of element blocks
C   --   NUMELB - IN - the number of elements per block
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   TIME - OUT - the time step time
C   --   WHOTIM - OUT - true iff whole (versus history) time step
C   --   VARHI - OUT - the history variables for the time step (if OPTION)
C   --   VARGL - OUT - the global variables for the time step (if OPTION)
C   --   VARNP - OUT - the nodal variables for the time step (if OPTION)
C   --   VAREL - OUT - the element variables for the time step (if OPTION)
C   --   * - return statement if error encountered, including end-of-file;
C   --      message is printed
C   --
C   --Database must be positioned in front of time step upon entry;
C   --upon exit positioned after time step.

      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER ISTEP
      INTEGER NVARHI, NVARGL, NVARNP, NUMNP, NVAREL, NELBLK
      INTEGER NUMELB(*)
c      LOGICAL ISEVOK(*)
      integer ISEVOK(*)
      REAL TIME
      LOGICAL WHOTIM
      REAL VARHI(*)
      REAL VARGL(*)
C     --NOTE: VARNP and VAREL are passed into DBIST? as doubly-dimensioned
      REAL VARNP(*)
      REAL VAREL(*)

      CHARACTER*80 ERRMSG

C   --Read step time

      READ (NDB, END=220, ERR=160, IOSTAT=IERR) TIME, HISTFL
      WHOTIM = (HISTFL .EQ. 0.0)

C   --Read history variables

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'H') .GT. 0)) THEN
         READ (NDB, END=170, ERR=170, IOSTAT=IERR)
     &      (VARHI(IVAR), IVAR=1,NVARHI)
      ELSE
         READ (NDB, END=170, ERR=170, IOSTAT=IERR)
      END IF

      IF (WHOTIM) THEN

C      --Read global variables

         IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'G') .GT. 0)) THEN
            READ (NDB, END=180, ERR=180, IOSTAT=IERR)
     &         (VARGL(IVAR), IVAR=1,NVARGL)
         ELSE
            READ (NDB, END=180, ERR=180, IOSTAT=IERR)
         END IF

C      --Read nodal variables

         IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'N') .GT. 0)) THEN
            CALL DBIST1 (NDB, NVARNP, NUMNP, VARNP, IVAR,
     &           *190)
         ELSE
            DO 110 IVAR = 1, NVARNP
               READ (NDB, END=190, ERR=190, IOSTAT=IERR)
  110       CONTINUE
         END IF

C      --Read element variables

         IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'E') .GT. 0)) THEN
            CALL DBIST2 (NDB, NVAREL,  NELBLK,
     $           MAX(NELBLK,1), ISEVOK, VAREL, NUMELB,
     &           IVAR, IELB, *200)
         ELSE
            DO 150 IELB = 1, NELBLK
               DO 140 IVAR = 1, NVAREL
                  IF ((ISEVOK( (IVAR-1)*NELBLK+IELB )) .ne. 0) THEN
                     READ (NDB, END=200, ERR=200, IOSTAT=IERR)
                  END IF
  140          CONTINUE
  150       CONTINUE
         END IF
      END IF

      RETURN

  160 CONTINUE
      WRITE (ERRMSG, '(A, I5)', IOSTAT=IDUM)
     &   'TIME for TIME STEP', ISTEP
      GOTO 210
  170 CONTINUE
      WRITE (ERRMSG, '(A, I5)', IOSTAT=IDUM)
     &   'HISTORY VARIABLES for TIME STEP', ISTEP
      GOTO 210
  180 CONTINUE
      WRITE (ERRMSG, '(A, I5)', IOSTAT=IDUM)
     &   'GLOBAL VARIABLES for TIME STEP', ISTEP
      GOTO 210
  190 CONTINUE
      WRITE (ERRMSG, '(A, I5, A, I5)', IOSTAT=IDUM)
     &   'NODAL VARIABLE', IVAR, ' for TIME STEP', ISTEP
      GOTO 210
  200 CONTINUE
      WRITE (ERRMSG, '(A, I5, A, I5, A, I5)', IOSTAT=IDUM)
     &   'ELEMENT VARIABLE', IVAR, ' of BLOCK', IELB,
     &   ' for TIME STEP', ISTEP
      GOTO 210
  210 CONTINUE
      CALL DBERR (IERR, ERRMSG)
  220 CONTINUE
      RETURN 1
      END
