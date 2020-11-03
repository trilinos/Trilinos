C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBITIM (NDB, OPTION, EXODUS,
     &   NVARNP, NELBLK, NVAREL, ISEVOK,
     &   NSTEPS, NSTEPW, A, KTIMES, KWHOLE, *)
C=======================================================================

C   --*** DBITIM *** (EXOLIB) Read database time step times
C   --   Written by Amy Gilkey - revised 11/11/87
C   --
C   --DBITIM reads all time steps of the database, storing the time for
C   --each step in dynamic memory.  It also sets the number of time steps.
C   --If there is a dynamic memory problem, the routine returns immediately
C   --without issuing an error message.  If there is a read error, a
C   --warning is printed and the step with the error is ignored.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   OPTION - IN - ' ' to not store, '*' to store all, else store options:
C   --      'W' to store whole times only
C   --      'H' to store history times only
C   --   EXODUS - IN - true iff this is an EXODUS file
C   --   NVARNP - IN - the number of nodal variables
C   --   NELBLK - IN - the number of element blocks
C   --   NVAREL - IN - the number of element variables
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   NSTEPS - OUT - the number of database time steps
C   --   NSTEPW - OUT - the number of whole (versus history) database time steps
C   --   A - IN/OUT - the dynamic memory base array
C   --   KTIMES - OUT - the dynamic memory index of TIMES (if OPTION)
C   --   KWHOLE - OUT - the dynamic memory index of WHOTIM (if OPTION);
C   --      WHOTIM(i) is true iff TIMES(i) is a whole (versus history) time step
C   --   * - return statement if error encountered, message is printed
C   --
C   --Database must be positioned in front of first time step upon entry;
C   --upon exit positioned at end of file.

      PARAMETER (MEMAMT=128)

      INTEGER NDB
      CHARACTER*(*) OPTION
      LOGICAL EXODUS
      INTEGER NVARNP, NELBLK, NVAREL
C     --NOTE: ISEVOK is doubly-dimensioned array (NELBLK,NVAREL)
      LOGICAL ISEVOK(*)
      INTEGER NSTEPS, NSTEPW
      REAL A(*)
      INTEGER KTIMES
      INTEGER KWHOLE

      LOGICAL STOWHO, STOHIS

      LOGICAL WHOTIM

      NSTEPS = 0
      NSTEPW = 0

      STOWHO = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'W') .GT. 0))
      STOHIS = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'H') .GT. 0))

      IF (STOWHO .OR. STOHIS) THEN
         NSTOR = 0

C      --Reserve memory for step times

         MAXMEM = MEMAMT
         CALL MDRSRV ('TIMES', KTIMES, MAXMEM)
         IF (STOWHO .AND. STOHIS) THEN
            CALL MDRSRV ('WHOTIM', KWHOLE, MAXMEM)
         END IF
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 150
      END IF

      IF (.NOT. EXODUS) GOTO 140

      WRITE (*, 10000) 'Please wait while the step times are read in'

C   --Set NRECST = the number of database records in a whole time step

      NRECST = 1 + 1 + 1 + NVARNP
      DO 110 IELB = 1, NELBLK
         DO 100 I = 1, NVAREL
            IF (ISEVOK( (I-1)*NELBLK+IELB) ) THEN
               NRECST = NRECST + 1
            END IF
  100    CONTINUE
  110 CONTINUE

  120 CONTINUE
      IF (.TRUE.) THEN

C      --Read step time

         READ (NDB, END=140, ERR=140, IOSTAT=IERR) TIME, HISTFL
         WHOTIM = (HISTFL .EQ. 0.0)

C      --Scan by variables

         IF (WHOTIM) THEN
            N = NRECST
         ELSE
            N = 2
         END IF
         DO 130 I = 2, N
            READ (NDB, END=160, ERR=160, IOSTAT=IERR)
  130    CONTINUE

         IF ((STOWHO .AND. WHOTIM) .OR.
     &      (STOHIS .AND. (.NOT. WHOTIM))) THEN

C         --Reserve more space for times if needed

            IF (NSTOR .GE. MAXMEM) THEN
               MAXMEM = MAXMEM + MEMAMT
               CALL MDLONG ('TIMES', KTIMES, MAXMEM)
               IF (STOWHO .AND. STOHIS) THEN
                  CALL MDLONG ('WHOTIM', KWHOLE, MAXMEM)
               END IF
               CALL MDSTAT (NERR, MEM)
               IF (NERR .GT. 0) GOTO 150
            END IF

C         --Store step time after time step read correctly

            A(KTIMES+NSTOR) = TIME
            IF (STOWHO .AND. STOHIS) THEN
               CALL CPYLOG (1, WHOTIM, A(KWHOLE+NSTOR))
            END IF

            NSTOR = NSTOR + 1
         END IF

         NSTEPS = NSTEPS + 1
         IF (WHOTIM) NSTEPW = NSTEPW + 1

         GOTO 120
      END IF

  140 CONTINUE

C   --Shorten space reserved for step times

      IF (STOWHO .OR. STOHIS) THEN
         CALL MDLONG ('TIMES', KTIMES, NSTOR)
         IF (STOWHO .AND. STOHIS) THEN
            CALL MDLONG ('WHOTIM', KWHOLE, NSTOR)
         END IF
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 150
      END IF

  150 CONTINUE
      RETURN

  160 CONTINUE
      CALL PRTERR ('WARNING', 'End-of-file during time steps')
      RETURN
10000  FORMAT (/, 1X, A)
      END
