C Copyright (c) 2007-2017 National Technology & Engineering Solutions of
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
C
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
C

C=======================================================================
*DECK,INTRPE
      SUBROUTINE INTRP3 (CNTRA,CNTRB,IELPT,SOLEB,SOLEA,SOLGRA,IDBLK,
     &                   ITT,iblk,TIMES,ISTP,IST,INSUB,ICOMPL,
     &                   XB,YB,ZB,ICONB,DUME)
C
C     ******************************************************************
C
C     SUBROUTINE TO CONTROL INTERPOLATION OF ELEMENT VARIABLES
C     FROM MESH-A TO MESH-B FOR SCHEME 3, A ELEMENT CENTROID BASED
C     INTERPOLATION SCHEME. PHYSICAL CONSTRAINTS ARE APPLIED TO
C     THE INTERPOLATED RESULTS AND THEN THEY ARE WRITTEN TO MESH-C
C
C     Called by MAPVAR
C
C     ******************************************************************
C
C CNTRA   REAL  Centroidal coords for Mesh-A
C CNTRB   REAL  Centroidal coords for Mesh-B
C IELPT   INT   The element in Mesh-A within which the point
C               (node in Mesh-B) is found
C STRPT   REAL  The isoparametric coords of point in IELPT element
C SOLEB   REAL  Element variables for Mesh-B
C SOLEA   REAL  Element variables for Mesh-A
C SOLGRA  REAL  Gradient in element variables for Mesh-A
C IDBLK   INT   Element block ID for mesh-B
C ITT     INT   Truth table
C IM      INT   Element block being processed (not block ID)
C TIMES   REAL  Array of times - passed thru to PNF
C ISTP    INT   Time step
C IST     INT   Time step if multiple time steps are in use
C INSUB   INT   Entry into subroutine; 1-first time in; >1-second,etc
C ICOMPL  INT   Map completion; 0-incomplete; 1-complete
C
C     ******************************************************************
C
      include 'aexds1.blk'
      include 'aexds2.blk'
      include 'amesh.blk'
      include 'bmesh.blk'
      include 'ebbyeb.blk'
      include 'ex2tp.blk'
      include 'tapes.blk'
C
      DIMENSION CNTRB(NUMEBB,*), CNTRA(NUMEBA,*), IELPT(*)
      DIMENSION SOLEB(NUMEBB,*), SOLEA(NUMEBA,*)
      DIMENSION SOLGRA(NDIMA,NUMEBA,*)
      DIMENSION ITT(NVAREL,*), ICONB(NELNDB,*)
      DIMENSION XX(27), YY(27), ZZ(27), XB(*), YB(*), ZB(*)
      DIMENSION DUME(*)
C
C     ******************************************************************
C
        IROT = 0
        IROTF = 0
        DO 40 IVAR=1,NVAREL
          IF (ITT(IVAR,iblk) .EQ. 0)GO TO 40
C
C Initialize SOLEB if first time in subroutine for this element block
C After first time into subroutine
C retrieve SOLEB from storage in EXODUS
C
        IF (INSUB .EQ. 1) THEN
          CALL INIELT(SOLEB,IVAR,TIMES,ISTP,IDBLK,CENTER,DUME)
        ELSE
          CALL EXGEV(NTP4EX,IST,IVAR,IDBLK,NUMEBB,SOLEB(1,IVAR),IERR)
        END IF
C
C Loop on elements in recipient mesh
C
        DO 30 I=1,NUMEBB
          IF (IELPT(I) .NE. 0)THEN
C
C Distance in cartesian coordinates between mesh-A and mesh-B centroid
C
            XC = CNTRB(I,1) - CNTRA(IELPT(I),1)
            YC = CNTRB(I,2) - CNTRA(IELPT(I),2)
            ZC = 0.
            IF (NDIMB .EQ. 3)ZC = CNTRB(I,3) - CNTRA(IELPT(I),3)
C
C Evaluate interpolation
C
            SOLEB(I,IVAR) = SOLEA(IELPT(I),IVAR)
     &                    + SOLGRA(1,IELPT(I),IVAR) * XC
     &                    + SOLGRA(2,IELPT(I),IVAR) * YC
     &                    + SOLGRA(3,IELPT(I),IVAR) * ZC
          END IF
 30     CONTINUE
C
C If there is more searching to do (i.e. many blocks to one)
C use EXODUS as temporary storage
C don't bother to perform needed adjustments yet
C
        IF (ICOMPL .NE. 1)THEN
          CALL EXPEV(NTP4EX,IST,IVAR,IDBLK,NUMEBB,SOLEB(1,IVAR),IERR)
        ELSE
C
C write element vars out to EXODUS data base (now is convenient)
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C ELMASS is special
C
           IF (NAMVAR(nvargp+IVAR)(1:6) .EQ. 'ELMASS')THEN
C
C ELMASS was changed to nodal density prior to processing.
C need to go back from density to element mass now
C
C             NNODES=NNELM(ITYPE)
             NNODES = NELNDB
             IF (ITYPE .EQ. 6) NNODES = 4
             DO 100 IEL = 1, NUMEBB
               DO 105 I = 1, NNODES
                 XX(I) = XB(ICONB(I,IEL))
                 YY(I) = YB(ICONB(I,IEL))
                 IF (NDIMB .EQ. 3)THEN
                   ZZ(I) = ZB(ICONB(I,IEL))
                 ELSE
                   ZZ(I) = 0.
                 END IF
 105           CONTINUE
               CALL VOL(ITYPE,XX,YY,ZZ,VOLUME)
               SOLEB(IEL,IVAR) = SOLEB(IEL,IVAR) * VOLUME
 100         CONTINUE
           END IF
C
C****************************************************************
C apply constraints to variables here as applicable
C
C
C Plastic strain (EQPS) must greater than or equal to 0.
C
           IF (NAMVAR(nvargp+IVAR)(1:4) .EQ. 'EQPS')THEN
             DO 110 IEL = 1, NUMEBB
               IF (SOLEB(IEL,IVAR) .LT. 0.)THEN
                 SOLEB(IEL,IVAR) = 0.
               END IF
  110        CONTINUE
C
C Hourglass forces and bulk viscosity have no meaning other than on
C the mesh from which they originated, just set them to zero.
C
           ELSE IF (NAMVAR(nvargp+IVAR)(1:2) .EQ. 'HG' .OR.
     &          NAMVAR(nvargp+IVAR)(1:5) .EQ. 'BULKQ')THEN
             DO 120 IEL = 1, NUMEBB
               SOLEB(IEL,IVAR) = 0.
  120        CONTINUE
           ELSE IF (NAMVAR(NVARGP+IVAR)(1:7) .EQ. 'TEARING')THEN
             DO 130 IEL = 1, NUMEBB
               IF (SOLEB(IEL,IVAR) .LT. 0.)THEN
                 SOLEB(IEL,IVAR) = 0.
               END IF
  130        CONTINUE
           END IF
C        ELSE IF (NAMVAR(NVARGP+IVAR)(1:?) .EQ. ?????)THEN
C          DO ??? IEL = 1, NUMEBB
C            IF (SOLEB(IEL,IVAR) .??. ?.)THEN
C              SOLEB(IEL,IVAR) = ?.
C            END IF
C  ???     CONTINUE
c************************************************************************
c
c########################################################################
c the rotation tensor is special
c
c just store pointers to the rotation tensor components for later
c processing. do nothing here
c
           IF (NAMVAR(NVARGP+IVAR)(1:8) .EQ. 'COSTHETA')THEN
             ICOS = IVAR
             IROT = IROT + 1
             GO TO 10
           ELSE IF (NAMVAR(NVARGP+IVAR)(1:8) .EQ. 'SINTHETA')THEN
             ISIN = IVAR
             IROT = IROT + 1
             GO TO 10
           ELSE IF (NAMVAR(NVARGP+IVAR)(1:3) .EQ. 'R11')THEN
             IR11 = IVAR
             IROT = IROT + 1
             GO TO 10
           ELSE IF (NAMVAR(NVARGP+IVAR)(1:3) .EQ. 'R21')THEN
             IR21 = IVAR
             IROT = IROT + 1
             GO TO 10
           ELSE IF (NAMVAR(NVARGP+IVAR)(1:3) .EQ. 'R31')THEN
             IR31 = IVAR
             IROT = IROT + 1
             GO TO 10
           ELSE IF (NAMVAR(NVARGP+IVAR)(1:3) .EQ. 'R12')THEN
             IR12 = IVAR
             IROT = IROT + 1
             GO TO 10
           ELSE IF (NAMVAR(NVARGP+IVAR)(1:3) .EQ. 'R22')THEN
             IR22 = IVAR
             IROT = IROT + 1
             GO TO 10
           ELSE IF (NAMVAR(NVARGP+IVAR)(1:3) .EQ. 'R32')THEN
             IR32 = IVAR
             IROT = IROT + 1
             GO TO 10
           ELSE IF (NAMVAR(NVARGP+IVAR)(1:3) .EQ. 'R13')THEN
             IR13 = IVAR
             IROT = IROT + 1
             GO TO 10
           ELSE IF (NAMVAR(NVARGP+IVAR)(1:3) .EQ. 'R23')THEN
             IR23 = IVAR
             IROT = IROT + 1
             GO TO 10
           ELSE IF (NAMVAR(NVARGP+IVAR)(1:3) .EQ. 'R33')THEN
             IR33 = IVAR
             IROT = IROT + 1
             GO TO 10
           END IF
c########################################################################
C
C write element variables
C
           CALL EXPEV(NTP4EX,IST,IVAR,IDBLK,NUMEBB,SOLEB(1,IVAR),IERR)
   10      CONTINUE
c
c########################################################################
c now fix-up rotations - rotation matrix must have mag=1
c
c some simple error checking
c
           IF (NDIMB .EQ. 2)THEN
             IF (IROT .EQ. 2 .AND.
     &         ICOS .NE. 0 .AND. ISIN .NE. 0) THEN
              DO 230 IEL = 1, NUMEBB
               RMAG = SQRT(SOLEB(IEL,ICOS)*SOLEB(IEL,ICOS)
     &                   + SOLEB(IEL,ISIN)*SOLEB(IEL,ISIN))
               SOLEB(IEL,ICOS) = SOLEB(IEL,ICOS) / RMAG
               SOLEB(IEL,ISIN) = SOLEB(IEL,ISIN) / RMAG
  230         CONTINUE
              CALL EXPEV(NTP4EX,IST,ICOS,IDBLK,NUMEBB,SOLEB(1,ICOS),
     &                   IERR)
              CALL EXPEV(NTP4EX,IST,ISIN,IDBLK,NUMEBB,SOLEB(1,ISIN),
     &                   IERR)
              IROTF = 1
             END IF
C
C
            ELSE IF (IROT .EQ. 9 .AND. IROTF .EQ. 0 .AND. IR11 .NE. 0
     &        .AND. IR21 .NE. 0 .AND. IR31 .NE. 0 .AND. IR12 .NE. 0
     &        .AND. IR22 .NE. 0 .AND. IR32 .NE. 0 .AND. IR13 .NE. 0
     &        .AND. IR23 .NE. 0 .AND. IR33 .NE. 0)THEN
C
C compute magnitude of matrix
C
            DO 280 IEL = I, NUMEBB
             RMAG=SQRT(SOLEB(IEL,IR11)*SOLEB(IEL,IR22)*SOLEB(IEL,IR33)
     &           +    SOLEB(IEL,IR21)*SOLEB(IEL,IR32)*SOLEB(IEL,IR13)
     &           +    SOLEB(IEL,IR31)*SOLEB(IEL,IR12)*SOLEB(IEL,IR23)
     &           -    SOLEB(IEL,IR11)*SOLEB(IEL,IR23)*SOLEB(IEL,IR32)
     &           -    SOLEB(IEL,IR12)*SOLEB(IEL,IR21)*SOLEB(IEL,IR33)
     &           -    SOLEB(IEL,IR13)*SOLEB(IEL,IR22)*SOLEB(IEL,IR31))

             SOLEB(IEL,IR11) = SOLEB(IEL,IR11) / RMAG
             SOLEB(IEL,IR21) = SOLEB(IEL,IR21) / RMAG
             SOLEB(IEL,IR31) = SOLEB(IEL,IR31) / RMAG
             SOLEB(IEL,IR12) = SOLEB(IEL,IR12) / RMAG
             SOLEB(IEL,IR22) = SOLEB(IEL,IR22) / RMAG
             SOLEB(IEL,IR32) = SOLEB(IEL,IR32) / RMAG
             SOLEB(IEL,IR13) = SOLEB(IEL,IR13) / RMAG
             SOLEB(IEL,IR23) = SOLEB(IEL,IR23) / RMAG
             SOLEB(IEL,IR33) = SOLEB(IEL,IR33) / RMAG
c
  280       CONTINUE
            CALL EXPEV(NTP4EX,IST,IR11,IDBLK,NUMEBB,SOLEB(1,IR11),
     &                 IERR)
            CALL EXPEV(NTP4EX,IST,IR21,IDBLK,NUMEBB,SOLEB(1,IR21),
     &                 IERR)
            CALL EXPEV(NTP4EX,IST,IR31,IDBLK,NUMEBB,SOLEB(1,IR31),
     &                 IERR)
            CALL EXPEV(NTP4EX,IST,IR12,IDBLK,NUMEBB,SOLEB(1,IR12),
     &                 IERR)
            CALL EXPEV(NTP4EX,IST,IR22,IDBLK,NUMEBB,SOLEB(1,IR22),
     &                 IERR)
            CALL EXPEV(NTP4EX,IST,IR32,IDBLK,NUMEBB,SOLEB(1,IR32),
     &                 IERR)
            CALL EXPEV(NTP4EX,IST,IR13,IDBLK,NUMEBB,SOLEB(1,IR13),
     &                 IERR)
            CALL EXPEV(NTP4EX,IST,IR23,IDBLK,NUMEBB,SOLEB(1,IR23),
     &                 IERR)
            CALL EXPEV(NTP4EX,IST,IR33,IDBLK,NUMEBB,SOLEB(1,IR33),
     &                 IERR)
            IROTF = 1
          END IF
        END IF
 40     CONTINUE
        IF (IROTF .NE. 1) THEN
          CALL ERROR('INTRPE',
     &      'ROTATION MATRIX NORMALLY REQUIRED FOR RESTART',
     &      'DIMENSION ',NDIMB,
     &      'NUMBER OF ROTATION MATRIX COMPONENTS FOUND',IROT,
     &      'THIS IS ONLY A WARNING',' ',0)
        END IF
c########################################################################
C
      RETURN
      END
