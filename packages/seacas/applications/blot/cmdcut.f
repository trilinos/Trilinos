C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CMDCUT (VERB, INLINE, IFLD, INTYP, CFIELD,
     &                   RFIELD, A, *)
C=======================================================================

C   --*** CMDCUT *** (MESH) Process cut commands
C   --   Written by Amy Gilkey - revised 03/11/88
C   --
C   --Parameters:
C   --   VERB    - I/O - the verbs for the SHOW command
C   --   INLINE  - I/O- the parsed input line for the log file
C   --   IFLD,   - I/O - the free-field reader index and fields
C   --   INTYP,  - I/O - the free-field reader index and fields
C   --   CFIELD, - I/O - the free-field reader index and fields
C   --   RFIELD  - I/O - the free-field reader index and fields
C   --   A       - IN  - the dynamic memory base array
C   --
C   --Common Variables:
C   --   Uses IS3DIM, NUMNPF of /D3NUMS/
C   --   Uses DFAC of /DEFORM/
C   --   Uses UNMESH of /MSHLIM/
C   --   Sets NEWCUT, ISCUT, CUTPT, CUTNRM of /CUTOPT/
C   --   Uses PKRMAT, PKRCEN of /PICK/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      include 'dbnums.blk'
      include 'd3nums.blk'
      include 'deform.blk'
      include 'mshlim.blk'
      include 'cutopt.blk'
      include 'pick.blk'

      CHARACTER*(*) VERB
      CHARACTER*(*) INLINE
      INTEGER     INTYP(*)
      CHARACTER*(*) CFIELD(*)
      REAL        RFIELD(*)
      DIMENSION A(*)

      LOGICAL ISON
      REAL RNUM(6)
      REAL CUTPLA(3,3)
      LOGICAL FFMATC
      LOGICAL LDUM1, LDUM2

C PUT VERB IN OUTPUT STRING

      CALL FFADDC (VERB, INLINE)

C CHECK THAT WE ARE ID 3D

      IF (.NOT. IS3DIM) THEN
         CALL PRTERR ('CMDERR', 'Command allowed in 3D only')
         GOTO 100
      END IF

C SET NEWCUT FLAG

      IF (ISCUT) NEWCUT = .TRUE.
      ISCUT = .FALSE.

C "CUT OFF" COMMAND

      IF (FFMATC (IFLD, INTYP, CFIELD, 'OFF', 3)) THEN
         CALL FFADDC ('OFF', INLINE)
         NEWCUT = .TRUE.
         ISCUT = .FALSE.

C ISSUE WARNING

      ELSE
         IF (DFAC .NE. 0.0) THEN
            CALL PRTERR ('CMDWARN',
     &           'Cut is performed on undeformed mesh')
         END IF

C "CUT SCREEN" COMMAND

         IF (FFMATC (IFLD, INTYP, CFIELD, 'SCREEN', 1)) THEN
            ISON = .TRUE.
C            --Pick two points forming a line with a third point forming
C            --a plane along the displayed Z axis

            CALL PICK2D ('first plane point', ISON,
     &            .TRUE., IFLD, INTYP, RFIELD,
     &            RNUM(1), RNUM(2), *100)
            CALL PICK2D ('other plane point', ISON,
     &            .FALSE., IFLD, INTYP, RFIELD,
     &            RNUM(3), RNUM(4), *100)

            CALL UNROT (1, 1, PKRMAT, PKRCEN,
     &            RNUM(1), RNUM(2), UNMESH(KNEA),
     &            CUTPLA(1,1), CUTPLA(2,1), CUTPLA(3,1))
            CALL UNROT (1, 1, PKRMAT, PKRCEN,
     &            RNUM(3), RNUM(4), UNMESH(KNEA),
     &            CUTPLA(1,2), CUTPLA(2,2), CUTPLA(3,2))
            CALL UNROT (1, 1, PKRMAT, PKRCEN,
     &            RNUM(1), RNUM(2), UNMESH(KFAR),
     &            CUTPLA(1,3), CUTPLA(2,3), CUTPLA(3,3))

C         --Determine which side is being cut away

            CALL PICK2D ('point in cut mesh', ISON,
     &            .FALSE., IFLD, INTYP, RFIELD,
     &            RNUM(5), RNUM(6), *100)
            CALL UNROT (1, 1, PKRMAT, PKRCEN,
     &            RNUM(5), RNUM(6), UNMESH(KNEA),
     &            RNUM(1), RNUM(2), RNUM(3))

C ENTER POINTS IN THE OUTPUT STRING

            CALL FFADDR (CUTPLA(1,1), INLINE)
            CALL FFADDR (CUTPLA(2,1), INLINE)
            CALL FFADDR (CUTPLA(3,1), INLINE)
            CALL FFADDR (CUTPLA(1,2), INLINE)
            CALL FFADDR (CUTPLA(2,2), INLINE)
            CALL FFADDR (CUTPLA(3,2), INLINE)
            CALL FFADDR (CUTPLA(1,3), INLINE)
            CALL FFADDR (CUTPLA(2,3), INLINE)
            CALL FFADDR (CUTPLA(3,3), INLINE)
            CALL FFADDR (RNUM(1), INLINE)
            CALL FFADDR (RNUM(2), INLINE)
            CALL FFADDR (RNUM(3), INLINE)

C GET CUT POINT AND NORMAL FROM THREE POINTS AND POINT IN MESH

            CALL PTSNRM(CUTPLA, RNUM, CUTPT, CUTNRM, IERR)
            IF(IERR .NE. 0) GO TO 100

C "CUT NORM" COMMAND

         ELSE IF (FFMATC (IFLD, INTYP, CFIELD, 'NORM', 1)) THEN
            CALL FFADDC ('NORM', INLINE)
            ISON =  FFMATC (IFLD, INTYP, CFIELD, 'CURSOR', 1)

C GET POINT ON CUT SURFACE

            CALL QNPICK ('DISPLAYED', LDUM1, LDUM2,
     &            A, KXN, KYN, KZN, KHIDEN, KNPSUR)
            CALL PICK3D ('point on cut surface', ISON,
     &            NUMNPF, A(KXN), A(KYN), A(KZN), A(KHIDEN),
     &            .TRUE., IFLD, INTYP, RFIELD,
     &            CUTPT(1), CUTPT(2), CUTPT(3), *100)

C GET POINT FOR NORMAL

            CALL PICK3D ('point for normal direction', ISON,
     &            NUMNPF, A(KXN), A(KYN), A(KZN), A(KHIDEN),
     &            .TRUE., IFLD, INTYP, RFIELD,
     &            CUTNRM(1), CUTNRM(2), CUTNRM(3), *100)

C ENTER CUT POINT IN OUTPUT STRING

            CALL FFADDR (CUTPT(1), INLINE)
            CALL FFADDR (CUTPT(2), INLINE)
            CALL FFADDR (CUTPT(3), INLINE)

C IF IN CURSOR MODE, THEN SUBTRACT NORMAL FROM CUT POINT TO GET THE
C NORMAL DIRECTION

            IF (ISON) THEN
               CUTNRM(1) = CUTNRM(1) - CUTPT(1)
               CUTNRM(2) = CUTNRM(2) - CUTPT(2)
               CUTNRM(3) = CUTNRM(3) - CUTPT(3)
            END IF

C ENTER NORMAL POINT IN OUTPUT STRING

            CALL FFADDR (CUTNRM(1), INLINE)
            CALL FFADDR (CUTNRM(2), INLINE)
            CALL FFADDR (CUTNRM(3), INLINE)

C CALCULATE NORMAL AND NORMALIZE

            DIST = SQRT(CUTNRM(1)*CUTNRM(1) + CUTNRM(2)*CUTNRM(2)
     &                  + CUTNRM(3)*CUTNRM(3))
            IF (DIST .EQ. 0) GO TO 100
            CUTNRM(1) = CUTNRM(1)/DIST
            CUTNRM(2) = CUTNRM(2)/DIST
            CUTNRM(3) = CUTNRM(3)/DIST

C "CUT" COMMANDS

         ELSE
            ISON =  FFMATC (IFLD, INTYP, CFIELD, 'CURSOR', 1)

C           --Pick three points on the cutting plane

            CALL QNPICK ('DISPLAYED', LDUM1, LDUM2,
     &            A, KXN, KYN, KZN, KHIDEN, KNPSUR)
            CALL PICK3D ('first plane point', ISON,
     &            NUMNPF, A(KXN), A(KYN), A(KZN), A(KHIDEN),
     &            .TRUE., IFLD, INTYP, RFIELD,
     &            CUTPLA(1,1), CUTPLA(2,1), CUTPLA(3,1), *100)
            CALL PICK3D ('second plane point', ISON,
     &            NUMNPF, A(KXN), A(KYN), A(KZN), A(KHIDEN),
     &            .FALSE., IFLD, INTYP, RFIELD,
     &            CUTPLA(1,2), CUTPLA(2,2), CUTPLA(3,2), *100)
            CALL PICK3D ('third plane point', ISON,
     &            NUMNPF, A(KXN), A(KYN), A(KZN), A(KHIDEN),
     &            .FALSE., IFLD, INTYP, RFIELD,
     &            CUTPLA(1,3), CUTPLA(2,3), CUTPLA(3,3), *100)

C         --Determine which side is being cut away

            CALL PICK3D ('point in cut mesh', ISON,
     &            NUMNPF, A(KXN), A(KYN), A(KZN), A(KHIDEN),
     &            .FALSE., IFLD, INTYP, RFIELD,
     &            RNUM(1), RNUM(2), RNUM(3), *100)

C ENTER POINTS IN OUTPUT STRING

            CALL FFADDR (CUTPLA(1,1), INLINE)
            CALL FFADDR (CUTPLA(2,1), INLINE)
            CALL FFADDR (CUTPLA(3,1), INLINE)
            CALL FFADDR (CUTPLA(1,2), INLINE)
            CALL FFADDR (CUTPLA(2,2), INLINE)
            CALL FFADDR (CUTPLA(3,2), INLINE)
            CALL FFADDR (CUTPLA(1,3), INLINE)
            CALL FFADDR (CUTPLA(2,3), INLINE)
            CALL FFADDR (CUTPLA(3,3), INLINE)
            CALL FFADDR (RNUM(1), INLINE)
            CALL FFADDR (RNUM(2), INLINE)
            CALL FFADDR (RNUM(3), INLINE)

C GET CUT POINT AND NORMAL FROM THREE POINTS AND POINT IN MESH

            CALL PTSNRM(CUTPLA, RNUM, CUTPT, CUTNRM, IERR)
            IF(IERR .NE. 0) GO TO 100

         END IF

C SET CUTTING FLAGS

         NEWCUT = .TRUE.
         ISCUT = .TRUE.
      END IF

      RETURN

  100 CONTINUE
      RETURN 1
      END
