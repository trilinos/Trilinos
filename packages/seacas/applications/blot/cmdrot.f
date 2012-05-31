C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
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
C     * Neither the name of Sandia Corporation nor the names of its
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

C=======================================================================
      SUBROUTINE CMDROT (VERB, INLINE,
     &   IFLD, INTYP, CFIELD, IFIELD, RFIELD,
     &   DFRMAT, DFRCEN,
     &   NEWZM, A, *)
C=======================================================================

C   --*** CMDROT *** (MESH) Process rotation commands
C   --   Written by Amy Gilkey - revised 03/07/88
C   --
C   --Parameters:
C   --   VERB - IN/OUT - the verbs for the SHOW command
C   --   INLINE - IN/OUT - the parsed input line for the log file
C   --   IFLD, INTYP, CFIELD, IFIELD, RFIELD - IN/OUT - the free-field
C   --      reader index and fields
C   --   DFRMAT, DFRCEN - IN - the default rotation or rotation center
C   --   NEWZM - IN - true iff a new zoom window or scaling is set
C   --   A - IN - the dynamic memory base array
C   --
C   --Common Variables:
C   --   Uses IS3DIM, NNPSUR, NUMNPF of /D3NUMS/
C   --   Uses UNMESH, RDMESH of /MSHLIM/
C   --   Sets and uses NEWROT, ROTMAT, ROTCEN, EYE of /ROTOPT/
C   --   Uses PKMESH, PKRMAT, PKRCEN of /PICK/

      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      include 'params.blk'
      include 'dbnums.blk'
      include 'd3nums.blk'
      include 'mshlim.blk'
      include 'rotopt.blk'
      include 'pick.blk'

      CHARACTER*(*) VERB
      CHARACTER*(*) INLINE
      INTEGER     INTYP(*)
      CHARACTER*(*) CFIELD(*)
      INTEGER     IFIELD(*)
      REAL        RFIELD(*)
      LOGICAL NEWZM
      REAL DFRMAT(3,3), DFRCEN(3)
      DIMENSION A(*)

      CHARACTER*(MXSTLN) WORD
      LOGICAL ISON
      REAL RNUM(3)
      LOGICAL FFEXST, FFMATC, MATSTR
      LOGICAL LDUM1, LDUM2

      IF (VERB .EQ. 'ROTATE') THEN
         CALL FFADDC (VERB, INLINE)
         IF (.NOT. IS3DIM) THEN
            CALL PRTERR ('CMDERR', 'Command allowed in 3D only')
            GOTO 120
         END IF

         DEGANG = 4.0 * ATAN(1.0) / 180.0
         NEWROT = .TRUE.

  100    CONTINUE
         IF (FFEXST (IFLD, INTYP)) THEN
            CALL FFCHAR (IFLD, INTYP, CFIELD, ' ', WORD)
            IF (MATSTR (WORD, 'RESET', 1)) THEN
               CALL FFADDC ('RESET', INLINE)
               CALL CPYREA (3*3, DFRMAT, ROTMAT)
            ELSE IF ((WORD .EQ. 'X') .OR. (WORD .EQ. 'Y')
     &         .OR. (WORD .EQ. 'Z')) THEN
               CALL FFREAL (IFLD, INTYP, RFIELD,
     &            'angle of rotation', 0.0, DEG, *110)
               CALL FFADDC (WORD, INLINE)
               CALL FFADDR (DEG, INLINE)
               CALL ROTXYZ (WORD, DEG * DEGANG, ROTMAT)
            ELSE
               CALL PRTERR ('CMDERR',
     &            'Expected "X", "Y", "Z" or "RESET"')
               GOTO 110
            END IF
            GOTO 100
         END IF

  110    CONTINUE
         Z = UNMESH(KFAR) - UNMESH(KNEA)
         CALL UNROT (1, 1, ROTMAT, ROTCEN,
     &      0.0, 0.0, Z, EYE(1), EYE(2), EYE(3))

      ELSE IF (VERB .EQ. 'EYE') THEN
         CALL FFADDC (VERB, INLINE)
         IF (.NOT. IS3DIM) THEN
            CALL PRTERR ('CMDERR', 'Command allowed in 3D only')
            GOTO 120
         END IF

         ISON = FFMATC (IFLD, INTYP, CFIELD, 'CURSOR', 1)
         CALL QNPICK ('DISPLAYED', LDUM1, LDUM2,
     &      A, KXN, KYN, KZN, KHIDEN, KNPSUR)
         CALL PICK3D ('eye position', ISON,
     &      NUMNPF, A(KXN), A(KYN), A(KZN), A(KHIDEN),
     &      .TRUE., IFLD, INTYP, RFIELD,
     &      RNUM(1), RNUM(2), RNUM(3), *120)
         CALL FFADDR (RNUM(1), INLINE)
         CALL FFADDR (RNUM(2), INLINE)
         CALL FFADDR (RNUM(3), INLINE)

         CALL ROTEYE (RNUM, ROTCEN, ROTMAT, *120)

         EYE(1) = RNUM(1)
         EYE(2) = RNUM(2)
         EYE(3) = RNUM(3)

         NEWROT = .TRUE.

      ELSE IF (VERB .EQ. 'CENTER') THEN
         CALL FFADDC (VERB, INLINE)
         IF (.NOT. IS3DIM) THEN
            CALL PRTERR ('CMDERR', 'Command allowed in 3D only')
            GOTO 120
         END IF

         IF (FFMATC (IFLD, INTYP, CFIELD, 'RESET', 1)) THEN
            CALL FFADDC ('RESET', INLINE)
            CALL CPYREA (3, DFRCEN, ROTCEN)

         ELSE IF (FFMATC (IFLD, INTYP, CFIELD, 'ZOOM', 1)) THEN
            CALL FFADDC ('ZOOM', INLINE)
            IF ((MSCTYP .NE. 'MESH') .AND. (MSCTYP .NE. 'ZOOM')) THEN
               CALL PRTERR ('CMDERR', 'ZOOM window is not defined')
               GOTO 120
            END IF

            IF (NEWZM) THEN
               CALL QNPICK ('ORIGINAL', LDUM1, LDUM2,
     &            A, KXN, KYN, KZN, KHIDEN, KNPSUR)
               CALL ROTZM (RDMESH,
     &            NNPSUR, A(KNPSUR), A(KXN), A(KYN), A(KZN),
     &            .TRUE., ROTMAT, ROTCEN, RDMESH, ROTCEN, *120)
            ELSE
               CALL QNPICK ('DISPLAYED', LDUM1, LDUM2,
     &            A, KXN, KYN, KZN, KHIDEN, KNPSUR)
               CALL ROTZM (PKMESH,
     &            NNPSUR, A(KNPSUR), A(KXN), A(KYN), A(KZN),
     &            .FALSE., PKRMAT, PKRCEN, RDMESH, ROTCEN, *120)
            END IF

         ELSE
            ISON = FFMATC (IFLD, INTYP, CFIELD, 'CURSOR', 1)
            CALL QNPICK ('DISPLAYED', LDUM1, LDUM2,
     &         A, KXN, KYN, KZN, KHIDEN, KNPSUR)
            CALL PICK3D ('center of rotation', ISON,
     &         NUMNPF, A(KXN), A(KYN), A(KZN), A(KHIDEN),
     &         .TRUE., IFLD, INTYP, RFIELD,
     &         RNUM(1), RNUM(2), RNUM(3), *120)
            CALL FFADDR (RNUM(1), INLINE)
            CALL FFADDR (RNUM(2), INLINE)
            CALL FFADDR (RNUM(3), INLINE)
            ROTCEN(1) = RNUM(1)
            ROTCEN(2) = RNUM(2)
            ROTCEN(3) = RNUM(3)
         END IF

         NEWROT = .TRUE.
      END IF

      RETURN

  120 CONTINUE
      RETURN 1
      END
