C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C========================================================================
*DECK, ELTON0
      SUBROUTINE ELTON0(ICONA,NELTN,SOLEA,SOLENA,
     &                  IDBLK,XA,YA,ZA,ISTP,ITT,iblk)

C  *********************************************************************

C  Subroutine ELTON0 extracts nodal values of element variables by
C  looping over each element and summing the value of the variable
C  at that element to each node in the connectivity list for that
C  element. Then the nodal summation of element variables is divided
C  by the number of elements that contributed to that node (resulting
C  in a nodal average of the element value.) This is done for the old
C  mesh elements and nodes to facilitate interpolation.

C  Each element block must be processed independently in order to
C  avoid averaging element variables across material boundaries.
C  Note: the last set of DO loops acts over all nodes; to make sense
C        one element block must be completely processed before another
C        element block is sent into this subroutine.

C  Calls subroutine VOL

c  Called by MAPVAR

C  *********************************************************************

C   ICONA       mesh-A connectivity (1:nelnda,1:numeba)
C   NELTN       number of elements tied to each node (1:nodesa)
C   SOLEA       element variables (1:numeba,1:nvarel)
C   SOLENA      element variables at nodes (1:nodesa,1:nvarel)
C   IDBLK       current element block I.D.
C   XA,YA,ZA    coordinates
C   XX,YY,ZZ    vector of coordinates of nodes for an element
C   ISTP        current time step
C   ITT         truth table
C   iblk        element block being processed (not ID)

C  *********************************************************************

      include 'aexds1.blk'
      include 'aexds2.blk'
      include 'amesh.blk'
      include 'ebbyeb.blk'
      include 'ex2tp.blk'
      include 'tapes.blk'

      DIMENSION ICONA(NELNDA,*), NELTN(*)
      DIMENSION SOLEA(NUMEBA,*), SOLENA(NODESA,NVAREL), ITT(NVAREL,*)
      DIMENSION XA(*), YA(*), ZA(*), XX(27), YY(27), ZZ(27)

C  *********************************************************************

      IF (ITYPE .EQ. 4 .OR. ITYPE .EQ. 5)THEN
        CALL ERROR('ELTON0','ELEMENT TYPE',' ',ITYPE,
     &             'ELEMENT VARIABLE PROCESSING NOT YET IMPLEMENTED',
     &              0,' ',' ',1)
      END IF

      DO I = 1, NODESA
         NELTN(I) = 0
         DO J = 1, NVAREL
            SOLENA(I,J) = 0.
         end do
      end do

C      NNODES = NNELM(ITYPE)
      NNODES = NELNDA
      IF (ITYPE .EQ. 6) NNODES = 4
      DO NEL = 1, NUMEBA
        DO I = 1, NNODES

C  number of elements associated with each node - used for
C    computing an average later on

          NELTN(ICONA(I,NEL)) = NELTN(ICONA(I,NEL)) + 1
       end do
      end do

      DO IVAR = 1, NVAREL
        IF (ITT(IVAR,iblk) .EQ. 0)GO TO 40
        CALL EXGEV(NTP2EX,ISTP,IVAR,IDBLK,NUMEBA,SOLEA(1,IVAR),IERR)

        IF (NAMVAR(nvargp+IVAR)(1:6) .EQ. 'ELMASS') THEN

C replace element mass with nodal density for interpolation

          DO IEL = 1, NUMEBA
            DO I = 1, NNODES
              XX(I) = XA(ICONA(I,IEL))
              YY(I) = YA(ICONA(I,IEL))
              IF (NDIMA .EQ. 3)THEN
                ZZ(I) = ZA(ICONA(I,IEL))
              ELSE
                ZZ(I) = 0.
              END IF
           end do
            CALL VOL(ITYPE,XX,YY,ZZ,VOLUME)
            SOLEA(IEL,IVAR) = SOLEA(IEL,IVAR) / VOLUME
         end do
        END IF

C  accumulate element variables to nodes

        DO NEL = 1, NUMEBA
          DO I = 1, NNODES
            SOLENA(ICONA(I,NEL),IVAR) =
     &      SOLENA(ICONA(I,NEL),IVAR) + SOLEA(NEL,IVAR)
         end do
      end do

C  divide by number of elements contributing to each node (average)

        DO I = 1, NODESA
          IF(NELTN(I) .NE. 0)THEN
            SOLENA(I,IVAR) = SOLENA(I,IVAR) / dble(NELTN(I))
          END IF
       end do
   40 CONTINUE
      end do
      RETURN
      END

