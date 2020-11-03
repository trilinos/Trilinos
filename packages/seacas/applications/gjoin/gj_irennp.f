C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C     -*- Mode: fortran -*-
C=======================================================================

      SUBROUTINE IRENNP (A, NNPS1, NNPS2, IDNPS, NNNPS,
     &     RENNP, MATNS1, MATNS2, TOLER, CLOSE, MATMAT,
     $     XSCL, YSCL, ZSCL, XOFF, YOFF, ZOFF, NDIM, IEXPCT)
C=======================================================================

C     --*** IRENNP *** (GJOIN) Input node renumber information from user
C     --   Written by Amy Gilkey - revised 02/23/88
C     --
C     --IRENNP requests node renumbering information from the user.  This
C     --includes whether the nodes need to be combined and if nodal point
C     --set should be used in combining the nodes.
C     --Other transformations to the mesh database are also performed here
C     --including: mirroring and offsetting
C     --
C     --Parameters:
C     --   A - IN/OUT - the dynamic memory array
C     --   NNPS1 - IN - the number of nodal point sets in the first database
C     --   NNPS2 - IN - the number of nodal point sets in the second database
C     --   IDNPS - IN - the IDs for each set
C     --   NNNPS - IN - the number of nodes for each set
C     --   RENNP - OUT - true iff the nodes need to be combined
C     --   MATNS1 - OUT - the number of the first nodal point set
C     --   MATNS2 - OUT - the number of the second nodal point set
C     --   TOLER - OUT - tolerance to be used for equivalencing
C     --   CLOSE - OUT - true if match closest node on equivalence
C                        (This is now default, cannot be turned off)

      PARAMETER (MAXFLD = 10)

      include 'gj_filnum.blk'
      include 'gj_xyzrot.blk'
      DIMENSION A(*)
      INTEGER IDNPS(*), NNNPS(*)
      LOGICAL RENNP, CLOSE, MATMAT

      CHARACTER*8 REPLY
      INTEGER     INTYP(MAXFLD+1)
      CHARACTER*8 CFIELD(MAXFLD), WORD
      INTEGER     IFIELD(MAXFLD)
      REAL        RFIELD(MAXFLD)
      CHARACTER*132 INPUT

      LOGICAL ONSET, MATSTR, FFEXST

      ROT3D   = .FALSE.
      CALL INIREA (3,   0.0, ROTCEN)
      CALL INIREA (3*3, 0.0, ROTMAT)
      DO 100 I = 1, 3
         ROTMAT(I,I) = 1.0
  100 continue
      IEXPCT  = 0
      CLOSE = .TRUE.
      MATMAT  = .FALSE.
      ONSET   = .FALSE.
      if (nnps1 + nnps2 .gt. 0) then
         CALL MDRSRV ('ISTA', KISTA, NNPS1+NNPS2)
         CALL MDRSRV ('ISCR', KISCR, NNPS1+NNPS2)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) THEN
            CALL MEMERR
            GOTO 240
         END IF
         CALL INITIN (A(KISTA), NNPS1+NNPS2, 0)
      end if
  110 CONTINUE
      WRITE (*, *)
      CALL GETINP (0,0,'Equivalence or Command (Enter HELP for info)? ',
     *  INPUT, IOSTAT)

C ... Handle commented lines
      IF (INPUT(1:1) .EQ. '$') goto 110

      IDCONT = 0
      CALL FFISTR(INPUT, MAXFLD, IDCONT,
     *  NUMFLD, INTYP, CFIELD, IFIELD, RFIELD)
      INTYP(MIN(MAXFLD,NUMFLD)+1) = -999

      CALL OUTLOG (KLOG, NUMFLD, INTYP, CFIELD, IFIELD, RFIELD)

C     -- Valid Commands:
C     NO
C     EQUIVALENCE
C     COMBINE
C     YES
C     LIST
C     MIRROR
C     OFFSET
C     SHIFT

      REPLY = CFIELD(1)
      IF (REPLY(1:1) .EQ. ' ') REPLY = 'NO'
      CALL EXUPCS (REPLY)

      IF (MATSTR(REPLY, 'HELP', 1)) THEN
         write (*, 10000)
10000    FORMAT (/,' Valid Commands:',/,4X,
     $        'YES|NO',/,4X,
     $        'EQUIVALENCE [nset1] [nset2] [toler]',
     $        ' [MATERIAL] ',/,4X,
     $        'COMBINE     [nset1] [nset2] [toler]',
     $        ' [MATERIAL]',/,4X,
     $        'EQUIVALENCE|COMBINE END|EXIT|NO',/,4X,
     $        'MIRROR X|Y|Z|ALL  (sets scale to -1)',/,4X,
     $        'OFFSET X|Y|Z|ALL|ADD|RESET offset ...',/,4X,
     $        'SHIFT  X|Y|Z|ALL|ADD|RESET offset ...',/,4X,
     $        'SCALE  X|Y|Z|ALL|MULTIPLY|RESET scale  ...',/,4X,
     $        'REVOLVE X|Y|Z|RESET angle ...',/,4x,
     $        'REVCEN xcen, ycen, zcen ',/,4x,
     $        'LIST (lists nodesets)',/,4X,
     $        'EXPECT num_matches',//,4X,
     $        'NOTE: Coords of 2nd database = scale*old_coord+offset',/)
         go to 110
      ELSE IF (MATSTR(REPLY, 'EQUIVALENCE', 3) .OR.
     $        MATSTR(REPLY, 'COMBINE', 3)) THEN
         IF (MATSTR(CFIELD(2), 'END', 3) .OR.
     $        MATSTR(CFIELD(2), 'EXIT',3) .OR.
     $        MATSTR(CFIELD(2), 'NO', 1)) THEN
            RENNP = .FALSE.
         ELSE
            RENNP = .TRUE.
         END IF
      ELSE IF (MATSTR(REPLY, 'YES', 1)) THEN
         RENNP = .TRUE.
      ELSE IF (MATSTR(REPLY, 'NO', 1) .or. MATSTR(REPLY, 'END',3)) THEN
         RENNP = .FALSE.
      else if (matstr(reply, 'LIST', 1)) then
         CALL PRTNPS (A(KISTA), NNPS1, NNPS2, IDNPS, NNNPS, A(KISCR))
         go to 110
      else if (matstr(reply, 'MIRROR', 1)) then
         ifld = 2
  120    CONTINUE
         IF (FFEXST (IFLD, INTYP)) THEN
            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (MATSTR (WORD, 'RESET', 1)) THEN
               XSCL = 1.
               YSCL = 1.
               ZSCL = 1.
            ELSE IF (WORD .EQ. 'X') THEN
               XSCL = -1.
            ELSE IF (WORD .EQ. 'Y') THEN
               YSCL = -1.
            ELSE IF (WORD .EQ. 'Z') THEN
               IF (NDIM .EQ. 3) THEN
                  ZSCL = -1.
               ELSE
                  CALL PRTERR ('CMDERR',
     *                 'Z not allowed for 2D database')
               END IF
            ELSE IF (WORD .EQ. 'ALL') THEN
               XSCL = -1.
               YSCL = -1.
               ZSCL = -1.
            ELSE
               IF (NDIM .EQ. 3) CALL PRTERR ('CMDERR',
     &              'Expected "X", "Y", "Z" or "RESET"')
               IF (NDIM .EQ. 2) CALL PRTERR ('CMDERR',
     &              'Expected "X", "Y", or "RESET"')
               GOTO 130
            END IF
            GOTO 120
         END IF

  130    CONTINUE
         go to 110
      ELSE IF (matstr(reply, 'REVOLVE', 4)) THEN
         ifld = 2
         DEGANG = ATAN2(0.0, -1.0) / 180.0

  140    CONTINUE
         IF (FFEXST (IFLD, INTYP)) THEN
            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (MATSTR (WORD, 'RESET', 1)) THEN
               ROT3D = .FALSE.
               CALL INIREA (3*3, 0.0, ROTMAT)
               DO 150 I = 1, 3
                  ROTMAT(I,I) = 1.0
  150          CONTINUE
            ELSE IF (NDIM .EQ. 3 .AND.
     *              ((WORD .EQ. 'X') .OR. (WORD .EQ. 'Y')
     &              .OR. (WORD .EQ. 'Z')) ) THEN
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &              'angle of rotation', 0.0, DEG, *160)
               ROT3D = .TRUE.
               CALL ROTXYZ (WORD, DEG * DEGANG, ROTMAT)
            ELSE IF (NDIM .EQ. 2 .AND. (WORD .EQ. 'Z') ) THEN
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &              'angle of rotation', 0.0, DEG, *160)
               ROT3D = .TRUE.
               CALL ROTXYZ (WORD, DEG * DEGANG, ROTMAT)
            ELSE
               IF (NDIM .EQ. 3) CALL PRTERR ('CMDERR',
     &              'Expected "X", "Y", "Z" or "RESET"')
               IF (NDIM .EQ. 2) CALL PRTERR ('CMDERR',
     &              'Expected "Z" or "RESET"')
               GOTO 160
            END IF
            GOTO 140
         END IF

  160    CONTINUE
         go to 110
      ELSE IF (matstr(reply, 'REVCEN', 4)) THEN
         ifld = 2
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &        'X revolution center', 0.0, ROTCEN(1), *170)
         CALL FFREAL (IFLD, INTYP, RFIELD,
     &        'Y revolution center', 0.0, ROTCEN(2), *170)
         IF (NDIM .EQ. 3) CALL FFREAL (IFLD, INTYP, RFIELD,
     &        'Z revolution center', 0.0, ROTCEN(3), *170)
         IF (NDIM .EQ. 2) ROTCEN(3) = 0.0

  170    CONTINUE
         go to 110
      else IF (matstr(reply, 'OFFSET', 1) .OR.
     $        matstr(reply, 'SHIFT',  2)) THEN
         ifld = 2
C     ... Originally offset just asked for NDIM values.  It was changed
C     to go by axis type (OFFSET Y 1.0). We do maintain compatibility
C     and check for both types of input.

         RMULT = 0.0
C         IF (INTYP(IFLD) .EQ. 0) THEN
C        Check if character field has 'X','Y', or 'Z'
         IF ((INTYP(IFLD) .EQ. 0) .OR. ((CFIELD(IFLD) .GE. 'A')
     &      .AND. (CFIELD(IFLD) .LE. 'Z'))) THEN
  180       CONTINUE
            IF (FFEXST (IFLD, INTYP)) THEN
               CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
               IF (MATSTR (WORD, 'RESET', 1)) THEN
                  XOFF = 0.0
                  YOFF = 0.0
                  ZOFF = 0.0
                  RMULT = 0.0
               ELSE IF (MATSTR (WORD, 'ADD', 2)) THEN
C     ... Set for cumulative offsets/shifts
                  RMULT = 1.0
               ELSE IF (MATSTR (WORD, 'ALL', 2)) THEN
                  CALL FFREAL (IFLD, INTYP, RFIELD,
     &                 'model offset', 0.0, TOFFS, *190)
                  XOFF = RMULT * XOFF + TOFFS
                  YOFF = RMULT * YOFF + TOFFS
                  ZOFF = RMULT * ZOFF + TOFFS
               ELSE IF (WORD .EQ. 'X') THEN
                  CALL FFREAL (IFLD, INTYP, RFIELD,
     &                 'X coordinate offset', 0.0, TOFFS, *190)
                  XOFF = RMULT * XOFF + TOFFS
               ELSE IF (WORD .EQ. 'Y') THEN
                  CALL FFREAL (IFLD, INTYP, RFIELD,
     &                 'Y coordinate offset', 0.0, TOFFS, *190)
                  YOFF = RMULT * YOFF + TOFFS
               ELSE IF (WORD .EQ. 'Z') THEN
                  IF (NDIM .EQ. 3) THEN
                     CALL FFREAL (IFLD, INTYP, RFIELD,
     &                    'Z coordinate offset', 0.0, TOFFS, *190)
                     ZOFF = RMULT * ZOFF + TOFFS
                  ELSE
                     CALL PRTERR ('CMDERR',
     *                    'Z allowed for 3D database only')
                  END IF
               ELSE
                  IF (NDIM .EQ. 3) CALL PRTERR ('CMDERR',
     &                 'Expected "X", "Y", "Z", "ALL", "ADD",'//
     $                 'or "RESET"')
                  IF (NDIM .EQ. 2) CALL PRTERR ('CMDERR',
     &                 'Expected "X", "Y", "ALL", "ADD", or "RESET"')
                  GOTO 190
               END IF
               GOTO 180
            END IF
         ELSE
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &           'X coordinate offset', XOFF, XOFF, *190)
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &           'Y coordinate offset', YOFF, YOFF, *190)
            IF (NDIM .EQ. 3) CALL FFREAL (IFLD, INTYP, RFIELD,
     &           'Z coordinate offset', ZOFF, ZOFF, *190)
         END IF
  190    CONTINUE
         go to 110
C     ------------------------------ SCALE -----------------------------
      else IF (matstr(reply, 'SCALE', 2)) THEN
         ifld = 2

         IF (INTYP(IFLD) .EQ. 0) THEN
  200       CONTINUE
            IF (FFEXST (IFLD, INTYP)) THEN
               CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
               IF (MATSTR (WORD, 'RESET', 1)) THEN
                  XSCL = 1.0
                  YSCL = 1.0
                  ZSCL = 1.0
               ELSE IF (MATSTR (WORD, 'ALL', 2)) THEN
                  CALL FFREAL (IFLD, INTYP, RFIELD,
     &                 'model scale', 1.0, TSCL, *210)
                  XSCL = TSCL
                  YSCL = TSCL
                  ZSCL = TSCL
               ELSE IF (WORD .EQ. 'X') THEN
                  CALL FFREAL (IFLD, INTYP, RFIELD,
     &                 'X coordinate scale', 1.0, XSCL, *210)
               ELSE IF (WORD .EQ. 'Y') THEN
                  CALL FFREAL (IFLD, INTYP, RFIELD,
     &                 'Y coordinate scale', 1.0, YSCL, *210)
               ELSE IF (WORD .EQ. 'Z') THEN
                  IF (NDIM .EQ. 3) THEN
                     CALL FFREAL (IFLD, INTYP, RFIELD,
     &                    'Z coordinate scale', 1.0, ZSCL, *210)
                  ELSE
                     CALL PRTERR ('CMDERR',
     *                    'Z allowed for 3D database only')
                  END IF
               ELSE
                  IF (NDIM .EQ. 3) CALL PRTERR ('CMDERR',
     &                 'Expected "X", "Y", "Z", "ALL",'//
     $                 'or "RESET"')
                  IF (NDIM .EQ. 2) CALL PRTERR ('CMDERR',
     &                 'Expected "X", "Y", "ALL", or "RESET"')
                  GOTO 210
               END IF
               GOTO 200
            END IF
         ELSE
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &           'X coordinate scale', XSCL, XSCL, *210)
            CALL FFREAL (IFLD, INTYP, RFIELD,
     &           'Y coordinate scale', YSCL, YSCL, *210)
            IF (NDIM .EQ. 3) CALL FFREAL (IFLD, INTYP, RFIELD,
     &           'Z coordinate scale', ZSCL, ZSCL, *210)
         END IF
  210    CONTINUE
         go to 110
C     ------------------------------ EXPECT -----------------------------
      else IF (matstr(reply, 'EXPECT', 3)) THEN
         ifld = 2
         call ffintg (ifld, intyp, ifield,
     $        'number of expected matches', 0, IEXPCT, *215)
 215     CONTINUE
         go to 110
      ELSE
         CALL PRTERR ('CMDERR',
     $        'Unrecognized response, Please Reenter')
         CALL PRTERR ('CMDSPEC',
     $        'Valid Response: Yes, No, COMbine, EQUivalence')
         GO TO 110
      END IF

      IF (IOSTAT .LT. 0) GOTO 240

      MATNS1 = 0
      MATNS2 = 0
      TOLER  = -999.0

C     ... Parse Command lines:
C     EQUIVALENCE [NS1 NS2] TOLER [CLOSEST|MATERIAL]

      IF (RENNP .AND. (NNPS1 .GT. 0) .AND. (NNPS2 .GT. 0)) THEN
         IF (NUMFLD .EQ. 2) THEN

C     ... 'EQUIV' TOLER
            TOLER   = RFIELD(2)

         ELSE IF (NUMFLD .EQ. 3 .AND. INTYP(3) .NE. 0) THEN

C     ... 'EQUIV' NS1 NS1
            ONSET   = .TRUE.

         ELSE IF (NUMFLD .EQ. 3 .AND. INTYP(3) .EQ. 0) THEN

C     ... 'EQUIV' TOLER 'CLOSEST' | 'MATERIAL'
            TOLER = RFIELD(2)
            IF (MATSTR(CFIELD(3),'CLOSEST',1)) THEN
               CLOSE = .TRUE.
            ELSE IF (MATSTR(CFIELD(3),'MATERIAL',1)) THEN
               MATMAT  = .TRUE.
            ELSE
               CALL PRTERR ('CMDWARN',
     $              'Unrecognized EQUIV Option: '//CFIELD(3))
            END IF

         ELSE IF (NUMFLD .EQ. 4 .and.
     $           intyp(3) .eq. 0 .and. intyp(4) .eq. 0) THEN

C     ... 'EQUIV' TOLER 'CLOSEST' 'MATERIAL'
            onset   = .FALSE.
            TOLER = RFIELD(2)
            IF (MATSTR(CFIELD(3),'CLOSEST',1))  CLOSE = .TRUE.
            IF (MATSTR(CFIELD(3),'MATERIAL',1)) MATMAT  = .TRUE.
            IF (MATSTR(CFIELD(4),'CLOSEST',1))  CLOSE = .TRUE.
            IF (MATSTR(CFIELD(4),'MATERIAL',1)) MATMAT  = .TRUE.

         ELSE IF (NUMFLD .EQ. 4 .and.
     $           intyp(3) .ne. 0 .and. intyp(4) .ne. 0) THEN
C     ... 'EQUIV' NS1 NS2 TOLER
            ONSET   = .TRUE.
            TOLER   = RFIELD(4)

         ELSE IF (NUMFLD .GE. 5) THEN

C     ... 'EQUIV' NS1 NS2 TOLER 'CLOSEST'|'MATERIAL'
            ONSET   = .TRUE.
            TOLER   = RFIELD(4)
            IF (MATSTR(CFIELD(5),'CLOSEST',1))  CLOSE = .TRUE.
            if (MATSTR(CFIELD(5),'MATERIAL',1))  MATMAT = .TRUE.
            if (numfld .eq. 6) then
               IF (MATSTR(CFIELD(6),'CLOSEST',1))  CLOSE = .TRUE.
               if (MATSTR(CFIELD(6),'MATERIAL',1))  MATMAT = .TRUE.
            end if
         ELSE IF (NUMFLD .EQ. 1) THEN
            CALL FREFLD (0, 0,
     *           'Should a nodal point set match be done? ',
     *           MAXFLD, IOSTAT, NUMFLD, INTYP, CFIELD, IFIELD, RFIELD)
            CALL OUTLOG (KLOG, NUMFLD, INTYP, CFIELD, IFIELD, RFIELD)

            REPLY = CFIELD(1)
            CALL EXUPCS (REPLY)
            ONSET = (REPLY(1:1) .EQ. 'Y')
            IF (IOSTAT .LT. 0) GOTO 240
         else
            call prterr ('PROGRAM',
     $           'Unrecognized command line format, try again.')
            write (*,*)
     $           'EQUIV [nset1 nset2] toler [material]'
         END IF

         IF (ONSET) THEN
            CALL PRTNPS (A(KISTA), NNPS1, NNPS2, IDNPS, NNNPS, A(KISCR))

  220       CONTINUE
            IF (NUMFLD .GE. 3) THEN
               IFLD = 2
               INTYP(MIN(MAXFLD,NUMFLD)+1) = -999
               NUMFLD = 1
            ELSE
               WRITE (*, *)
               CALL FREFLD (0, 0,
     &              'Enter set ID of first set, second set> ', MAXFLD,
     &              IOSTAT, NUMFLD, INTYP, CFIELD, IFIELD, RFIELD)
               IF (IOSTAT .LT. 0) GOTO 240
               IF (NUMFLD .EQ. 0) GOTO 220
               CALL OUTLOG (KLOG, NUMFLD, INTYP, CFIELD, IFIELD, RFIELD)

               INTYP(MIN(MAXFLD,NUMFLD)+1) = -999
               IFLD = 1
            END IF
            CALL FFINTG (IFLD, INTYP, IFIELD,
     &           'first set ID', 0, ID1, *220)
            CALL FFINTG (IFLD, INTYP, IFIELD,
     &           'second set ID', 0, ID2, *220)
            IX1 = LOCINT (ID1, NNPS1, IDNPS)
            IX2 = LOCINT (ID2, NNPS2, IDNPS(NNPS1+1))
            IF ((IX1 .LE. 0) .AND. (IX2 .LE. 0)) THEN
               CALL PRTERR ('CMDERR', 'Invalid set ID for both sets')
            ELSE IF (IX1 .LE. 0) THEN
               CALL PRTERR ('CMDERR', 'Invalid set ID for first set')
            ELSE IF (IX2 .LE. 0) THEN
               CALL PRTERR ('CMDERR', 'Invalid set ID for second set')
            END IF
            IF ((IX1 .LE. 0) .OR. (IX2 .LE. 0)) GOTO 110
            IX2 = NNPS1 + IX2
            MATNS1 = IX1
            MATNS2 = IX2

         END IF
      ELSE IF (NUMFLD .GE. 2) THEN
         TOLER = RFIELD(2)
         DO 230 IFLD = 3, NUMFLD
            IF (MATSTR(CFIELD(IFLD),'CLOSEST',1)) THEN
               CLOSE = .TRUE.
            ELSE IF (MATSTR(CFIELD(IFLD),'MATERIAL',1)) THEN
               MATMAT  = .TRUE.
            ELSE
               CALL PRTERR ('CMDWARN',
     $              'Unrecognized EQUIV Option: '//CFIELD(IFLD))
            END IF
  230    CONTINUE
      END IF
  240 CONTINUE
      if (nnps1 + nnps2 .gt. 0) then
         CALL MDDEL ('ISTA')
         CALL MDDEL ('ISCR')
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) CALL MEMERR
      end if

      RETURN
      END
