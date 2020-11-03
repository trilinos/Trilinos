C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
*DECK,INTRPN
      SUBROUTINE INTRPN(ICONA,SOLNA,IELPT,STRPT,
     &                  SOLNB,NDLSTB,XB,YB,ZB,
     &                  IDBLK,TIMES,INSUB,SN)

C     ******************************************************************

C     SUBROUTINE TO CONTROL INTERPOLATION OF NODAL RESULTS FROM MESH-A
C     TO MESH-B
C     INTERPOLATED SOLUTION IS WRITTEN TO MESH-C EXODUS FILE

C     Calls subroutine SHAPEF, ININOD

C     Called by MAPVAR

C     ******************************************************************

C ICONA   INT   Connectivity of donor Mesh (1:nelnda,1:numeba)
C SOLNA   REAL  Nodal variables  for donor mesh
C IELPT   INT   The element in donor mesh within which the point
C               (node in recipient mesh) is found
C STRPT   REAL  The isoparametric coords of point in IELPT element
C SOLNB   REAL  Nodal variables for recipient mesh
C NDLSTB  INT   List of recipient mesh nodes in current element block
C SOLN    REAL  SOLNA vector local to each donor mesh element
C TIMES   REAL  Array of times - just passed through
C INSUB   INT   Number of times into subroutine for this element block
C               1-first time; >1-second,etc time;
C               used to control mapping for many element blocks to one

C     ******************************************************************

      include 'amesh.blk'
      include 'bmesh.blk'
      include 'aexds1.blk'
      include 'contrl.blk'
      include 'ebbyeb.blk'
      include 'ex2tp.blk'
      include 'steps.blk'
      include 'tapes.blk'
      include 'varnpt.blk'

      DIMENSION XB(*), YB(*), ZB(*), TIMES(*)
      DIMENSION ICONA(NELNDA,*), SOLNA(NODESA,NVARNP)
      DIMENSION SOLNB(NODESB,NVARNP), NDLSTB(*)
      DIMENSION IELPT(*), STRPT(3,NODESB)
      DIMENSION SOLN(27), SN(*)

C     ******************************************************************

C set up time steps

      IF (ISTEP .EQ. -1)THEN
        NTM = NTIMES
      ELSE
        NTM = 1
      END IF

      DO 5 IST = 1, NTM
        IF (ISTEP .EQ. -1)THEN
          ISTP = IST
        ELSE
          ISTP = ISTEP
        END IF

C Start interpolation

        DO 10 IVAR = 1, NVARNP

C For IDEF = 2 do mesh annealing (used by GOMA); write out all
C displacements as zero.

          IF (IDEF .EQ. 2 .AND. (IVAR .EQ. IXDIS .OR. IVAR .EQ. IYDIS
     &        .OR. IVAR .EQ. IZDIS))GO TO 10

C If first time into INTRPN for this element block, initialize
C else you are mapping many to one and retrieve partially mapped
C results from temporary storage in EXODUS

          IF (INSUB .EQ. 1)THEN
            CALL ININOD(SOLNB,IVAR,TIMES,ISTP,IDBLK,NDLSTB,XB,YB,ZB,
     &                  SN)
          ELSE
            CALL EXGNV(NTP4EX,IST,IVAR,NODESB,SOLNB(1,IVAR),IERR)
          END IF

C Get nodal results on donor mesh

          CALL EXGNV(NTP2EX,ISTP,IVAR,NODESA,SOLNA(1,IVAR),IERR)

C Loop on nodes in recipient mesh

          DO 30 I = 1,NUMNDB
            NEL = IELPT(I)
            IF (NEL .NE. 0) THEN

C Set parameters for element in donor mesh

              S = STRPT(1,I)
              T = STRPT(2,I)
              R = 0.
              IF (NDIMB.EQ.3) R = STRPT(3,I)
C              NNODES = NNELM(ITYPE)
              NNODES = NELNDA
              IF (ITYPE .EQ. 6) NNODES = 4
              DO 20 J = 1,NNODES
                INODE = ICONA(J,NEL)
                SOLN(J) = SOLNA(INODE,IVAR)
 20           CONTINUE

C Shape function

              CALL SHAPEF(ITYPE,S,T,R,SOLN,BVALUE)
              SOLNB(NDLSTB(I),IVAR) = BVALUE
            END IF
 30       CONTINUE

C Save results, it doesn't matter if they are preliminary or final

          CALL EXPNV(NTP4EX,IST,IVAR,NODESB,SOLNB(1,IVAR),IERR)

 10     CONTINUE
  5   CONTINUE

      RETURN
      END
