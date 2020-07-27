C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
*DECK,INTRPE
      SUBROUTINE INTRPE (ICONA,SOLENA,IELPT,STRPT,
     &                   SOLEB,IDBLK,XB,YB,ZB,
     &                   ICONB,ITT,IBLK, TIMES,
     &                   CENTER,ISTP,IST,INSUB,ICOMPL,DUME)

C     ******************************************************************

C     SUBROUTINE TO CONTROL INTERPOLATION OF ELEMENT TRANSFORMED
C     INTO NODAL RESULTS FROM MESH-A TO MESH-B
C     INTERPOLATED SOLUTION IS PASSED OUT TO BE RETRANSFORMED
C     INTO ELEMENT RESULTS AND THEN WRITTEN TO MESH-C EXODUS FILE

C     Calls subroutines SHAPEF

C     Called by MAPVAR

C     ******************************************************************

C ICONA   INT   Connectivity of Mesh-A (1:nelnda,1:numeba)
C SOLENA  REAL  Element variables at nodes for Mesh-A
C IELPT   INT   The element in Mesh-A within which the point
C               (node in Mesh-B) is found
C STRPT   REAL  The isoparametric coords of point in IELPT element
C SOLEB   REAL  Element variables for Mesh-B
C SOLN    REAL  SOLENA vector local to each Mesh-A element
C TIMES   REAL  Array of times - passed through to PNF
C CENTER  REAL  Coordinates of element centroid - passed through to PNF
C ISTP    INT   Time step
C INSUB   INT   Entry into subroutine; 1-first time in; >1-second,etc
C ICOMPL  INT   Map completion; 0-incomplete; 1-complete

C     ******************************************************************

      include 'aexds1.blk'
      include 'aexds2.blk'
      include 'amesh.blk'
      include 'bmesh.blk'
      include 'ebbyeb.blk'
      include 'ex2tp.blk'
      include 'tapes.blk'

      DIMENSION TIMES(*), CENTER(NUMEBB,*)
      DIMENSION ICONA(NELNDA,*), SOLENA(NODESA,NVAREL)
      DIMENSION ITT(NVAREL,*)
      DIMENSION SOLEB(NUMEBB,*), XB(*), YB(*), ZB(*)
      DIMENSION IELPT(*), STRPT(3,NODESB), ICONB(NELNDB,*)
      DIMENSION SOLN(27), XX(27), YY(27), ZZ(27), DUME(*)

C     ******************************************************************

        IROT = 0
        IROTF = 0
        DO 40 IVAR=1,NVAREL
          IF (ITT(IVAR,IBLK) .EQ. 0)GO TO 40

C Initialize SOLEB if first time in subroutine for this element block
C After first time into subroutine
C retrieve SOLEB from storage in EXODUS

        IF (INSUB .EQ. 1) THEN
          CALL INIELT(SOLEB,IVAR,TIMES,ISTP,IDBLK,CENTER,DUME)
        ELSE
          CALL EXGEV(NTP4EX,IST,IVAR,IDBLK,NUMEBB,SOLEB(1,IVAR),IERR)
        END IF

C Loop on centroids in recipient mesh

        DO 30 I=1,NUMEBB
          IF (IELPT(I) .NE. 0)THEN

C Set parameters for element in donor mesh

            S=STRPT(1,I)
            T=STRPT(2,I)
            R=0.
            IF (NDIMB.EQ.3) R=STRPT(3,I)
C            NNODES=NNELM(ITYPE)
            NNODES = NELNDB
            IF (ITYPE .EQ. 6) NNODES = 4
            DO 20 J=1,NNODES
              INODE=ICONA(J,IELPT(I))
              SOLN(J)=SOLENA(INODE,IVAR)
 20         CONTINUE

C Shape function to evaluate interpolation

            CALL SHAPEF (ITYPE,S,T,R,SOLN,BVALUE)
            SOLEB(I,IVAR) = BVALUE
          END IF
 30     CONTINUE

C If there is more searching to do (i.e. many blocks to one)
C use EXODUS as temporary storage
C don't bother to perform needed adjustments yet

        IF (ICOMPL .NE. 1)THEN
          CALL EXPEV(NTP4EX,IST,IVAR,IDBLK,NUMEBB,SOLEB(1,IVAR),IERR)
        ELSE

C write element vars out to EXODUS data base (now is convenient)

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C ELMASS is special

           IF (NAMVAR(nvargp+IVAR)(1:6) .EQ. 'ELMASS')THEN

C ELMASS was changed to nodal density prior to processing.
C need to go back from density to element mass now

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

C****************************************************************
C apply constraints to variables here as applicable

C Plastic strain (EQPS) must greater than or equal to 0.

           IF (NAMVAR(nvargp+IVAR)(1:4) .EQ. 'EQPS')THEN
             DO 110 IEL = 1, NUMEBB
               IF (SOLEB(IEL,IVAR) .LT. 0.)THEN
                 SOLEB(IEL,IVAR) = 0.
               END IF
  110        CONTINUE

C Hourglass forces and bulk viscosity have no meaning other than on
C the mesh from which they originated, just set them to zero.

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

c########################################################################
c the rotation tensor is special

c just store pointers to the rotation tensor components for later
c processing. do nothing here

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

C write element variables

           CALL EXPEV(NTP4EX,IST,IVAR,IDBLK,NUMEBB,SOLEB(1,IVAR),IERR)
   10      CONTINUE

c########################################################################
c now fix-up rotations - rotation matrix must have mag=1

c some simple error checking

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

            ELSE IF (IROT .EQ. 9 .AND. IROTF .EQ. 0 .AND. IR11 .NE. 0
     &        .AND. IR21 .NE. 0 .AND. IR31 .NE. 0 .AND. IR12 .NE. 0
     &        .AND. IR22 .NE. 0 .AND. IR32 .NE. 0 .AND. IR13 .NE. 0
     &        .AND. IR23 .NE. 0 .AND. IR33 .NE. 0)THEN

C compute magnitude of matrix

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
c########################################################################

      RETURN
      END
