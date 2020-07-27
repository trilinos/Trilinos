C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBIV0 (NQAREC, NINFO)
C=======================================================================
C   --*** DBIV0 *** (EXOLIB) Initialize for DBIVAR
C   --
C   --DBIV0 and DBIV1 initialize for the DBIVAR routine.
C   --
C   --Parameters:
C   --   NQAREC - IN - the number of QA records
C   --   NINFO - IN - the number of information records
C   --
C   --Database state is maintained within DBIVIN and DBIVAR, and should
C   --not be moved between calls to these routines.

      INTEGER NQAREC, NINFO

      INTEGER NELBLK
      INTEGER NVARHI, NVARGL, NVARNP, NVAREL
c      LOGICAL ISEVOK(*)
      integer ISEVOK(*)

      LOGICAL REWDB

      INTEGER NDB
      INTEGER NUMVAR
      INTEGER IVAR(*)
      INTEGER ISTEP
      INTEGER LENVAR
      INTEGER NUMELB(*)
      REAL VAR(*)

      CHARACTER*80 ERRMSG
      CHARACTER TYP, T
      LOGICAL WHOTIM

      INTEGER NQASV, NINSV, NELBSV, NVELSV
      SAVE NQASV, NINSV, NELBSV, NVELSV
C      --NQASV, NINSV, NELBSV, NVELSV - the saved values for
C      --   NQAREC, NINFO, NELBLK, NVAREL

      INTEGER NREC0, NRECST
      SAVE NREC0, NRECST
C      --NREC0 - the number of database records before the time steps
C      --NRECST - the number of database records for each whole time step

      INTEGER IHVR0, IGVR0, INVR0, IEVR0
      SAVE IHVR0, IGVR0, INVR0, IEVR0
C      --IHVR0, IGVR0, INVR0, IEVR0 - the variable record number offset,
C      --   including the time record

      INTEGER NCSTEP, NCREC, NCEND
      SAVE NCSTEP, NCREC, NCEND
C      --NCSTEP - the current database step number; <0 if database should be
C      --   rewound
C      --NCREC - the current database record within time step
C      --NCEND - the number of database records in current time step

      DATA NQASV, NINSV, NELBSV / -999, -999, -999 /
      DATA NREC0, NRECST / -999, -999 /
      DATA IHVR0, IGVR0, INVR0, IEVR0 / -999, -999, -999, -999 /
      DATA NCSTEPL / -999 /

C   --Save the input parameters

      NQASV = NQAREC
      NINSV = NINFO

      GOTO 120

C=======================================================================
      ENTRY DBIV1 (NELBLK, NVARHI, NVARGL, NVARNP, NVAREL, ISEVOK)
C=======================================================================

C   --*** DBIV1 *** (EXOLIB) Initialize for DBIVAR
C   --   Written by Amy Gilkey - revised 11/04/87
C   --
C   --DBIV0 and DBIV1 initialize for the DBIVAR routine.
C   --
C   --Parameters:
C   --   NELBLK - IN - the number of element blocks
C   --   NVARHI - IN - the number of history variables
C   --   NVARGL - IN - the number of global variables
C   --   NVARNP - IN - the number of nodal variables
C   --   NVAREL - IN - the number of element variables
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --
C   --Database state is maintained within DBIVIN and DBIVAR, and should
C   --not be moved between calls to these routines.

C   --Save the input parameters

      NELBSV = NELBLK
      NVELSV = NVAREL

C   --Set NRECST = the number of records in a whole time step

      IF (NRECST .LT. 0) THEN
         NRECST = 1 + 1 + 1 + NVARNP
         DO 110 IELB = 1, NELBLK
            DO 100 I = 1, NVAREL
               IF (ISEVOK(0+(ielb-1)*NVAREL+i) .ne. 0) THEN
                  NRECST = NRECST + 1
               END IF
  100       CONTINUE
  110    CONTINUE
      END IF

C   --Set the variable record number offsets

      IF (IHVR0 .LT. 0) THEN
         IHVR0 = 1
         IGVR0 = IHVR0 + 1
         INVR0 = IGVR0 + 1
         IEVR0 = INVR0 + NVARNP
      END IF

C   --DBIV0 and DBIV1 converge at this point
  120 CONTINUE

C   --Set NREC0 = the number of records before the time steps

      IF (NREC0 .LT. 0) THEN
         IF ((NELBSV .GE. 0)
     &      .AND. (NQASV .GE. 0) .AND. (NINSV .GE. 0)) THEN
            NREC0 = 1 + 1 + 1 + 1 + NELBSV * (1 + 1 + 1) + 5 + 8
     &         + 1 + NQASV + 1 + NINSV + 1 + 1 + 1 + 1 + 1
         END IF
      END IF

      RETURN

C=======================================================================
      ENTRY DBIVIN (REWDB)
C=======================================================================

C   --*** DBIVIN *** (EXOLIB) Initialize for DBIVAR
C   --   Written by Amy Gilkey - revised 11/04/87
C   --
C   --DBIVIN initializes for the DBIVAR routine.
C   --
C   --Parameters:
C   --   REWDB - IN - if true, database should be rewound, else database is
C   --      at start of time steps
C   --
C   --Database state is maintained within DBIVIN and DBIVAR, and should
C   --not be moved between calls to these routines.

      IF (REWDB) THEN
         NCSTEP = -999
      ELSE
         NCSTEP = 1
         NCREC = 0
      END IF

      RETURN

C=======================================================================
      ENTRY DBIVAR (NDB, NUMVAR, IVAR, ISTEP, LENVAR, IELBLK,
     &   NELBLK, NUMELB, ISEVOK, VAR, *)
C=======================================================================

C   --*** DBIVAR *** (EXOLIB) Read variable
C   --   Written by Amy Gilkey - revised 08/16/88
C   --
C   --DBIVAR returns the values for the requested variables for the
C   --requested time step.  Several variables may be read, but all must
C   --be of the same type.  Element variables may be read for a single
C   --block (because of the way they are stored, this prevents rewinding).
C   --
C   --DBIVIN must be called before this routine is called and the database
C   --cannot be moved between calls to DBIVIN and DBIVAR.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   NUMVAR - IN - the number of variables to be read (for nodal and
C   --      element variables only)
C   --   IVAR - IN - the variable indices; no repetitions allowed
C   --   ISTEP - IN - the time step number
C   --   LENVAR - IN - the length of VAR
C   --   IELBLK - IN - the element block number to read, 0 for all
C   --      (for element variables only)
C   --   NELBLK - IN - the number of element blocks
C   --      (for element variables only)
C   --   NUMELB - IN - the number of elements per block
C   --      (for element variables only)
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(i,j)
C   --      (for element variables only)
C   --      NOTE: ISEVOK is indexed as singly-dimensioned array.
C   --   VAR - OUT - the variable values
C   --   * - return statement if error encountered, message is printed

C   --Find the variable type and ID

      CALL DBVTYP (IVAR(1), TYP, MINID)
      MAXID = MINID
      DO 130 I = 2, NUMVAR
         CALL DBVTYP (IVAR(I), T, ID)
         IF (TYP .NE. T) THEN
            CALL PRTERR ('PROGRAM',
     &         'Variables of different types requested')
            GOTO 320
         END IF
         MINID = MIN (MINID, ID)
         MAXID = MAX (MAXID, ID)
  130 CONTINUE

      IF (MINID .EQ. 0) TYP = 'T'
      IF (TYP .EQ. ' ') THEN
         CALL PRTERR ('PROGRAM', 'Invalid variable requested')
         GOTO 320
      END IF

C   --Find the starting variable record number

      IF (TYP .EQ. 'T') THEN
         IVREC = 1
      ELSE IF (TYP .EQ. 'H') THEN
         IVREC = IHVR0 + 1
      ELSE IF (TYP .EQ. 'G') THEN
         IVREC = IGVR0 + 1
      ELSE IF (TYP .EQ. 'N') THEN
         IVREC = INVR0 + MINID
      ELSE IF (TYP .EQ. 'E') THEN
         IF (IELBLK .GE. 1) THEN
            ISELB = IELBLK
            IEELB = IELBLK
         ELSE
            ISELB = 1
            IEELB = NELBLK
         END IF

C      --Check that some record must be read, otherwise return
         DO 150 IELB = ISELB, IEELB
            DO 140 I = 1, NUMVAR
               CALL DBVTYP (IVAR(I), T, ID)
               IF (ISEVOK(0+(IELB-1)*NVELSV+id) .ne. 0) GOTO 160
  140       CONTINUE
  150    CONTINUE
         RETURN
  160    CONTINUE

         IVREC = IEVR0 + 1
         DO 180 IELB = 1, ISELB-1
            DO 170 I = 1, NVELSV
               IF (ISEVOK(0+(IELB-1)*NVELSV+i) .ne. 0) THEN
                  IVREC = IVREC + 1
               END IF
  170       CONTINUE
  180    CONTINUE
         DO 190 ID = 1, MINID-1
            IF (ISEVOK(0+(ISELB-1)*NVELSV+ID) .ne. 0) THEN
               IVREC = IVREC + 1
            END IF
  190    CONTINUE
      END IF

C   --Rewind the database if past record (or state is unknown)

      IF ((NCSTEP .LT. 0) .OR. (ISTEP .LT. NCSTEP)
     &   .OR. ((ISTEP .EQ. NCSTEP) .AND. (IVREC .LE. NCREC))) THEN
         REWIND (NDB, ERR=310)

         NCSTEP = 0
         NCREC = 0
         DO 200 I = NCREC+1, NREC0
            READ (NDB, END=310, ERR=310, IOSTAT=IERR)
            NCREC = NCREC + 1
  200    CONTINUE
         NCSTEP = 1
         NCREC = 0
      END IF

C   --Scan rest of partially-read current time step (if not needed)

      IF ((NCSTEP .LT. ISTEP) .AND. (NCREC .GT. 0)) THEN
         DO 210 I = NCREC+1, NCEND
            READ (NDB, END=310, ERR=310, IOSTAT=IERR)
            NCREC = NCREC + 1
  210    CONTINUE
         NCSTEP = NCSTEP + 1
         NCREC = 0
      END IF

C   --Scan past preceding time steps

      DO 230 IS = NCSTEP, ISTEP-1
         NCREC = 0
         READ (NDB, END=310, ERR=310, IOSTAT=IERR) TIME, HISTFL
         NCREC = NCREC + 1
         WHOTIM = (HISTFL .EQ. 0.0)

         IF (WHOTIM) THEN
            NCEND = NRECST
         ELSE
            NCEND = 2
         END IF
         DO 220 I = NCREC+1, NCEND
            READ (NDB, END=310, ERR=310, IOSTAT=IERR)
            NCREC = NCREC + 1
  220    CONTINUE
         NCSTEP = NCSTEP + 1
         NCREC = 0
  230 CONTINUE

C   --Scan to needed record in this time step

      IF (NCREC .LE. 0) THEN
         READ (NDB, END=310, ERR=310, IOSTAT=IERR) TIME, HISTFL
         NCREC = 1
         WHOTIM = (HISTFL .EQ. 0.0)

         IF (WHOTIM) THEN
            NCEND = NRECST
         ELSE
            NCEND = 2
         END IF
      END IF

      IF (IVREC .GT. NCEND) THEN
         CALL PRTERR ('PROGRAM',
     &      'History variables only on this time step')
         GOTO 320
      END IF
      DO 240 I = NCREC+1, IVREC-1
         READ (NDB, END=310, ERR=310, IOSTAT=IERR)
         NCREC = NCREC + 1
  240 CONTINUE

C   --Read variable values

      IF (TYP .EQ. 'T') THEN
         VAR(1) = TIME

      ELSE IF ((TYP .EQ. 'H') .OR. (TYP .EQ. 'G')) THEN
         READ (NDB, END=310, ERR=310, IOSTAT=IERR)
     &      (VAR(I), I=1,LENVAR)
         NCREC = NCREC + 1

      ELSE IF (TYP .EQ. 'N') THEN
         CALL DBVIX (TYP, 1, IX1)
         DO 250 ID = MINID, MAXID
            IX = LOCINT (ID+IX1-1, NUMVAR, IVAR)
            IF (IX .GT. 0) THEN
               READ (NDB, END=310, ERR=310, IOSTAT=IERR)
     &            (VAR( (IX-1)*LENVAR+I ), I=1,LENVAR)
               NCREC = NCREC + 1
            END IF
  250    CONTINUE

      ELSE IF (TYP .EQ. 'E') THEN
         CALL DBVIX (TYP, 1, IX1)

         IEL0 = 0
         DO 260 IELB = 1, ISELB-1
            IEL0 = IEL0 + NUMELB(IELB)
  260    CONTINUE

         DO 300 IELB = ISELB, IEELB
            IF (IELB .GT. ISELB) THEN
               DO 270 ID = 1, MINID-1
                  IF (ISEVOK( 0+(IELB-1)*NVELSV+ID) .ne. 0) THEN
                     READ (NDB, END=310, ERR=310, IOSTAT=IERR)
                     NCREC = NCREC + 1
                  END IF
  270          CONTINUE
            END IF

            DO 280 ID = MINID, MAXID
               IX = LOCINT (ID+IX1-1, NUMVAR, IVAR)
               IF (IX .GT. 0) THEN
                  IF (ISEVOK( 0+(IELB-1)*NVELSV+id) .ne. 0) THEN
                     READ (NDB, END=310, ERR=310, IOSTAT=IERR)
     &                  (VAR( (IX-1)*LENVAR+IEL0+NE), NE=1,NUMELB(IELB))
                     NCREC = NCREC + 1
                  END IF
               ELSE
                  IF (ISEVOK( 0+(IELB-1)*NVELSV+ID) .ne. 0) then
                     read (ndb, end=310, err=310, iostat=ierr)
                     ncrec = ncrec + 1
                  end if
               end if
  280       continue

            if (ielb .lt. ieelb) then
               do 290 id = maxid+1, nvelsv
                  if (isevok( (id-1)*nelblk+ielb ) .ne. 0) then
                     read (ndb, end=310, err=310, iostat=ierr)
                     ncrec = ncrec + 1
                  end if
  290          continue
            end if
            iel0 = iel0 + numelb(ielb)
  300    continue
      end if

      if (ncrec .ge. ncend) then
         ncstep = ncstep + 1
         ncrec = 0
      end if

      return

  310 continue
      errmsg = 'time step variable'
      call dberr (ierr, errmsg)
      ncstep = -999
  320 continue
      return 1
      end
