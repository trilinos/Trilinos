C Copyright(C) 2009-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C     * Neither the name of NTESS nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

C $Log: cmdcut.f,v $
C Revision 1.3  2009/03/25 12:36:42  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  1996/06/21 16:07:03  caforsy
C Ran ftnchek and removed unused variables.  Reformat output for list
C var, list global, and list name.
C
C Revision 1.1  1994/04/07 19:55:57  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:48:14  gdsjaar
c Added RCS Id and Log to all files
c
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

C
C PUT VERB IN OUTPUT STRING
C
      CALL FFADDC (VERB, INLINE)
C
C CHECK THAT WE ARE ID 3D
C
      IF (.NOT. IS3DIM) THEN
         CALL PRTERR ('CMDERR', 'Command allowed in 3D only')
         GOTO 100
      END IF
C
C SET NEWCUT FLAG
C
      IF (ISCUT) NEWCUT = .TRUE.
      ISCUT = .FALSE.
C
C "CUT OFF" COMMAND
C
      IF (FFMATC (IFLD, INTYP, CFIELD, 'OFF', 3)) THEN
         CALL FFADDC ('OFF', INLINE)
         NEWCUT = .TRUE.
         ISCUT = .FALSE.
C
C ISSUE WARNING
C
      ELSE
         IF (DFAC .NE. 0.0) THEN
            CALL PRTERR ('CMDWARN',
     &           'Cut is performed on undeformed mesh')
         END IF

C
C "CUT SCREEN" COMMAND
C
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
C
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
C
C GET CUT POINT AND NORMAL FROM THREE POINTS AND POINT IN MESH
C
            CALL PTSNRM(CUTPLA, RNUM, CUTPT, CUTNRM, IERR)
            IF(IERR .NE. 0) GO TO 100
C
C "CUT NORM" COMMAND
C
         ELSE IF (FFMATC (IFLD, INTYP, CFIELD, 'NORM', 1)) THEN
            CALL FFADDC ('NORM', INLINE)
            ISON =  FFMATC (IFLD, INTYP, CFIELD, 'CURSOR', 1)
C
C GET POINT ON CUT SURFACE
C
            CALL QNPICK ('DISPLAYED', LDUM1, LDUM2,
     &            A, KXN, KYN, KZN, KHIDEN, KNPSUR)
            CALL PICK3D ('point on cut surface', ISON,
     &            NUMNPF, A(KXN), A(KYN), A(KZN), A(KHIDEN),
     &            .TRUE., IFLD, INTYP, RFIELD,
     &            CUTPT(1), CUTPT(2), CUTPT(3), *100)
C
C GET POINT FOR NORMAL
C
            CALL PICK3D ('point for normal direction', ISON,
     &            NUMNPF, A(KXN), A(KYN), A(KZN), A(KHIDEN),
     &            .TRUE., IFLD, INTYP, RFIELD,
     &            CUTNRM(1), CUTNRM(2), CUTNRM(3), *100)
C
C ENTER CUT POINT IN OUTPUT STRING
C
            CALL FFADDR (CUTPT(1), INLINE)
            CALL FFADDR (CUTPT(2), INLINE)
            CALL FFADDR (CUTPT(3), INLINE)
C
C IF IN CURSOR MODE, THEN SUBTRACT NORMAL FROM CUT POINT TO GET THE
C NORMAL DIRECTION
C
            IF (ISON) THEN
               CUTNRM(1) = CUTNRM(1) - CUTPT(1)
               CUTNRM(2) = CUTNRM(2) - CUTPT(2)
               CUTNRM(3) = CUTNRM(3) - CUTPT(3)
            END IF
C
C ENTER NORMAL POINT IN OUTPUT STRING
C
            CALL FFADDR (CUTNRM(1), INLINE)
            CALL FFADDR (CUTNRM(2), INLINE)
            CALL FFADDR (CUTNRM(3), INLINE)
C
C CALCULATE NORMAL AND NORMALIZE
C

            DIST = SQRT(CUTNRM(1)*CUTNRM(1) + CUTNRM(2)*CUTNRM(2)
     &                  + CUTNRM(3)*CUTNRM(3))
            IF (DIST .EQ. 0) GO TO 100
            CUTNRM(1) = CUTNRM(1)/DIST
            CUTNRM(2) = CUTNRM(2)/DIST
            CUTNRM(3) = CUTNRM(3)/DIST
C
C "CUT" COMMANDS
C
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
C
C ENTER POINTS IN OUTPUT STRING
C
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

C
C GET CUT POINT AND NORMAL FROM THREE POINTS AND POINT IN MESH
C
            CALL PTSNRM(CUTPLA, RNUM, CUTPT, CUTNRM, IERR)
            IF(IERR .NE. 0) GO TO 100

         END IF
C
C SET CUTTING FLAGS
C
         NEWCUT = .TRUE.
         ISCUT = .TRUE.
      END IF

      RETURN

  100 CONTINUE
      RETURN 1
      END
