C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=====================================================================
*     DECK,MKE
      SUBROUTINE MKE(NELND,NUMEB,NUMND,ICON,NDLST,ITYPE,
     &     VELX,VELY,VELZ,EMSS,RNMS,
     &     RMX,RMY,RMZ,RKE,PSQ,RJ2,
     &     SIGXX,SIGYY,SIGZZ,SIGXY,SIGYZ,SIGZX)
C
C     ****************************************************************
C
C     Compute momenta (x, y, and z) and kinetic energy
C     for an element block. Totals are stored in MKEI.
C     Called when check on accuracy of map is requested
C
C     Called by MAPVAR
C
C     ****************************************************************
C
C     NELND  INT   Number of nodes per element
C     NUMEB  INT   Number of elements in element block
C     NUMND  int   Number of nodes in element block
C     ICON   INT   Connectivity array (1:nelnd,1:numeb)
C     NDLST  INT   Array of nodes in element block (1:nneb)
C     ITYPE  INT   Element type
C     VELX   REAL  Array of velocities from data base (1:nnod)
C     VELY   REAL  Array of velocities from data base (1:nnod)
C     VELZ   REAL  Array of velocities from data base (1:nnod)
C     EMSS   REAL  Array of element masses from data base (1:numeb)
C     RNMS   REAL  Array of nodal masses (1:nnod)
C     RMX    REAL  X-momentum for this element block, this time step
C     RMY    REAL  Y-momentum for this element block, this time step
C     RMZ    REAL  Z-momentum for this element block, this time step
C     RKE    REAL  Kinetic energy for this element block, this time step
C     PSQ    REAL  Pressure squared * element mass for element block
C     RJ2    REAL  J2 * element mass for element block
C     SIGXX  REAL  Component of stress tensot
C     SIGYY  REAL  Component of stress tensot
C     SIGZZ  REAL  Component of stress tensot
C     SIGXY  REAL  Component of stress tensot
C     SIGYZ  REAL  Component of stress tensot
C     SIGZX  REAL  Component of stress tensot

C
C     ****************************************************************
C
      include 'amesh.blk'
C
      DIMENSION ICON(NELND,*),NDLST(*)
      DIMENSION VELX(*),VELY(*),VELZ(*),EMSS(*),RNMS(*)
      DIMENSION SIGXX(*),SIGYY(*),SIGZZ(*),SIGXY(*),SIGYZ(*),SIGZX(*)
C
C     ****************************************************************
C
C     zero nodal mass array
C
      DO I = 1, NUMEB
         DO J = 1, NELND
            RNMS(ICON(J,I)) = 0.
         end do
      end do
C
C     Translate element mass to nodal mass
C     First cut, come back and do better if necessary
C
      IF (ITYPE .EQ. 3 .OR. ITYPE .EQ. 4 .OR. ITYPE .EQ. 5)THEN
C
C     Treat all quads the same, use only four corner nodes
C     can fix if needed
C
         NEND = 4
         DO I = 1, NUMEB
            DO J = 1, NEND
               RNMS(ICON(J,I)) = RNMS(ICON(J,I)) + 0.25 * EMSS(I)
            end do
         end do
      ELSE IF (ITYPE .EQ. 10)THEN
C
C     8-node hex
C
         NEND = 8
         DO I = 1, NUMEB
            DO J = 1, NEND
               RNMS(ICON(J,I)) = RNMS(ICON(J,I)) + 0.125 * EMSS(I)
            end do
         end do
      ELSE IF (ITYPE .EQ. 13)THEN
C
C     4-node shell
C
         NEND = 4
         DO I = 1, NUMEB
            DO J = 1, NEND
               RNMS(ICON(J,I)) = RNMS(ICON(J,I)) + 0.25 * EMSS(I)
            end do
         end do
      ELSE IF (ITYPE .EQ. 6)THEN
C
C     treat all tets the same, use only four corner nodes
C     can fix if needed
C
         NEND = 4
         DO I = 1, NUMEB
            DO J = 1, NEND
               RNMS(ICON(J,I)) = RNMS(ICON(J,I)) + 0.25 * EMSS(I)
            end do
         end do
      ELSE
         CALL ERROR ('MKE','UNSUPPORTED ELEMENT TYPE','ITYPE',itype,
     &        ' ',0,' ',' ',1)
      END IF
C
C     compute momenta, KE, P-sqM, J2-sqM for element block
C     loop over all nodes in element block
C
      RMX = 0.
      RMY = 0.
      RMZ = 0.
      RKE = 0.
      PSQ = 0.
      RJ2 = 0.
      DO I = 1, NUMND
         INOD = NDLST(I)
         IF (NDIMA .EQ. 3)THEN
            VELSQ = VELX(INOD)*VELX(INOD) + VELY(INOD)*VELY(INOD)
     &           + VElZ(INOD)*VELZ(INOD)
            RMZ = RMZ + (VELZ(INOD) * RNMS(INOD))
         ELSE
            VELSQ = VELX(INOD)*VELX(INOD) + VELY(INOD)*VELY(INOD)
         END IF
         RMX = RMX + (VELX(INOD) * RNMS(INOD))
         RMY = RMY + (VELY(INOD) * RNMS(INOD))
         RKE = RKE + (0.5 * RNMS(INOD) * VELSQ)
      end do
C
      DO I = 1, NUMEB
         P = (SIGXX(I) + SIGYY(I) + SIGZZ(I)) / 3.
         PSQ = PSQ + P*P*EMSS(I)
         TAUP1 = SIGXX(I) - P
         TAUP2 = SIGYY(I) - P
         TAUP3 = SIGZZ(I) - P
         IF (NDIMA .EQ. 3)THEN
            RJ2 = RJ2 + ((((TAUP1**2 + TAUP2**2 + TAUP3**2)/6.)
     &           + SIGXY(I)**2 + SIGYZ(I)**2 + SIGZX(I)**2) * EMSS(I))
         ELSE
            RJ2 = RJ2 + ((((TAUP1**2 + TAUP2**2 + TAUP3**2)/6.)
     &           + SIGXY(I)**2) * EMSS(I))
         END IF
      end do
C
      RETURN
      END
