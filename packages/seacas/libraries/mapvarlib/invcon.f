C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C========================================================================
      SUBROUTINE INVCON(INVLN,MAXLN,INVCN,ICONA)
C
C************************************************************************
C
C     Subroutine INVC0N computes the inverse connectivity (elements connected
C     to a node).
C
c     Called by MAPVAR
C
C     Calls ERROR
C
C************************************************************************
C
C     INVLN  INT   The number of elements connected to a node (1:numnda)
C     MAXLN  INT   The maximum number of elements connected to any node
C     INVCN  INT   The inverse connectivity (1:maxln,1:numnda)
C     ICONA  INT   The connectivity array (1:nelnda,1:numela)
C
C************************************************************************
C
C
      include 'amesh.blk'
      include 'ebbyeb.blk'
C
      DIMENSION INVLN(*),INVCN(MAXLN,*),ICONA(nelnda,*)
C
C************************************************************************
C
      DO I = 1, NODESA
         INVLN(I) = 0
         DO J = 1, MAXLN
            INVCN(J,I) = 0
         end do
      end do

C
      NNODES = NELNDA
      IF (ITYPE .EQ. 6) NNODES = 4
      DO J = 1, NUMEBA
         DO I = 1, NNODES
            node = icona(i,j)
            IF (invln(node) .eq. 0 .or.
     *           INVCN(INVLN(node),node) .NE. J) THEN
               INVLN(node) = INVLN(node) + 1
               IF (INVLN(node) .GT. MAXLN)
     &              CALL ERROR('INVCON',' ',
     &              'TOO MANY ELEMENTS CONNECTED TO NODE',
     &              node,'INVCN ARRAY DIMENSIONED FOR NO MORE THAN',
     &              MAXLN,'RESET IN SUBROUTINE RDA2',' ',1)
               INVCN(INVLN(node),node) = J
            END IF
         end do
      end do
      RETURN
      END

