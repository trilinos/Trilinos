C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
*DECK, ELGRAD
      SUBROUTINE ELGRAD(CNTRA,XA,YA,ZA,SOLEA,SOLGRA,ICHKEL,IDBLK,
     &                  ICONA,INVLN,INVCN,MAXLN,ISTP,ITT,iblk)

C  *********************************************************************

C  Subroutine ELGRAD computes the coefficients of the gradient for
C  each variable within an element based on a constrained (to the
C  value of that element at its centroid) least squares fit to the
C  values in the surrounding elements. If not enough data is available
C  to perform the least squares is a direction, that coefficient is set
C  to zero.

C  Calls subroutines CNTR, VOL, FLGRAD, ERROR

C  Called by MAPVAR

C  *********************************************************************

C   CNTRA     a list of element centroid coordinates for all elements
C             in the current element block (1:numeba,1:ndima)
C   XA,YA,ZA  nodal coordinates
C   SOLEA     element variables (1:numeba,1:nvarel)
C   SOLGRA    element variable gradient (1:ndima,1:numeba,1:nvarel)
C   ICHKEL    check if element already in IELLST
C   IDBLK     current element block I.D.
C   ICONA     connectivity (1:nelnda,1:numeba)
C   INVLN     number of elements per node (1:numnda)
C   INVCN     inverse connectivity (1:maxln,1:numnda)
C   MAXLN     maximum number of elements connected to any node
C   ISTP      time step being processed
C   IELLST    local list of elements that share a node with element
C             currently being processed (1:100)
C   SHLNRM    vector normal to current shell element (1:3)
C   ITT       truth table
C   iblk      element block being processed (not ID)

C  *********************************************************************

      include 'aexds1.blk'
      include 'aexds2.blk'
      include 'amesh.blk'
      include 'ebbyeb.blk'
      include 'ex2tp.blk'
      include 'tapes.blk'

      DIMENSION CNTRA(NUMEBA,*), SOLEA(NUMEBA,*)
      DIMENSION SOLGRA(NDIMA,NUMEBA,*), ICHKEL(*)
      DIMENSION XX(27), YY(27), ZZ(27), IELLST(100), SHLNRM(3)
      DIMENSION XA(*), YA(*), ZA(*), ICONA(NELNDA,*)
      DIMENSION INVCN(MAXLN,*),INVLN(*), ITT(NVAREL,*)

C  *********************************************************************

      IF (ITYPE .EQ. 4 .OR. ITYPE .EQ. 5)THEN
        CALL ERROR('ELGRAD','ELEMENT TYPE',' ',ITYPE,
     &             'ELEMENT VARIABLE PROCESSING NOT YET IMPLEMENTED',
     &              0,' ',' ',1)
      END IF

C initialize ICHKEL array - only needs to be done once because
C checkmark used is specific to each element

      DO 10 I = 1, NUMEBA
        ICHKEL(I) = 0
 10   CONTINUE

C  load up CNTRA array - coordinates of mesh-A element centroids

C      NNODES = NNELM(ITYPE)
      NNODES = NELNDA
      IF (ITYPE .EQ. 6) NNODES = 4
      IF (NDIMA .EQ. 2) THEN
        DO 20 IEL = 1, NUMEBA
          DO 30 I = 1, NNODES
            INODE = ICONA(I,IEL)
            XX(I) = XA(INODE)
            YY(I) = YA(INODE)
            ZZ(I) = 0.
   30     CONTINUE
          CALL CNTR(ITYPE,XX,YY,ZZ,CNTRA(IEL,1),CNTRA(IEL,2),DUMMY)

   20   CONTINUE
      ELSE
        DO 40 IEL = 1, NUMEBA
          DO 50 I = 1, NNODES
            INODE = ICONA(I,IEL)
            XX(I) = XA(INODE)
            YY(I) = YA(INODE)
            ZZ(I) = ZA(INODE)
   50     CONTINUE
          CALL CNTR(ITYPE,XX,YY,ZZ,CNTRA(IEL,1),CNTRA(IEL,2),
     &              CNTRA(IEL,3))

   40   CONTINUE
      END IF

C  put element variables into SOLEA array

      DO 80 IVAR = 1, NVAREL
        IF (ITT(IVAR,iblk) .EQ. 0)GO TO 80
        CALL EXGEV(NTP2EX,ISTP,IVAR,IDBLK,NUMEBA,SOLEA(1,IVAR),IERR)

        IF (NAMVAR(nvargp+IVAR)(1:6) .EQ. 'ELMASS') THEN

C replace element mass with density

          DO 90 IEL = 1, NUMEBA
            DO 100 I = 1, NNODES
              XX(I) = XA(ICONA(I,IEL))
              YY(I) = YA(ICONA(I,IEL))
              IF (NDIMA .EQ. 3)THEN
                ZZ(I) = ZA(ICONA(I,IEL))
              ELSE
                ZZ(I) = 0.
              END IF
  100       CONTINUE
            CALL VOL(ITYPE,XX,YY,ZZ,VOLUME)

            SOLEA(IEL,IVAR) = SOLEA(IEL,IVAR) / VOLUME
   90     CONTINUE
        END IF
   80 CONTINUE

C start element centroid based gradient calculation

C get list of elements that share at least one node with current
C element being processed IELLST

      DO 110 IEL = 1, NUMEBA
        ICOUNT = 0
        ICHKEL(IEL) = IEL
        DO 120 INODE = 1, NNODES
          NOWNOD = ICONA(INODE,IEL)
          DO 130 IE = 1, INVLN(NOWNOD)
            NOWELT = INVCN(IE,NOWNOD)
            IF (IEL .NE. ICHKEL(NOWELT))THEN
              ICOUNT = ICOUNT + 1
              ICHKEL(NOWELT) = IEL
              IELLST(ICOUNT) = NOWELT
            END IF
 130      CONTINUE
 120    CONTINUE

C If shells element, compute element normal
C use cross product of vectors from mid-sides

        IF (ITYPE .EQ. 13)THEN
          XM1 = XA(ICONA(1,IEL)) - XA(ICONA(2,IEL))
          YM1 = YA(ICONA(1,IEL)) - YA(ICONA(2,IEL))
          ZM1 = ZA(ICONA(1,IEL)) - ZA(ICONA(2,IEL))
          XM2 = XA(ICONA(2,IEL)) - XA(ICONA(3,IEL))
          YM2 = YA(ICONA(2,IEL)) - YA(ICONA(3,IEL))
          ZM2 = ZA(ICONA(2,IEL)) - ZA(ICONA(3,IEL))
          XM3 = XA(ICONA(3,IEL)) - XA(ICONA(4,IEL))
          YM3 = YA(ICONA(3,IEL)) - YA(ICONA(4,IEL))
          ZM3 = ZA(ICONA(3,IEL)) - ZA(ICONA(4,IEL))
          XM4 = XA(ICONA(4,IEL)) - XA(ICONA(1,IEL))
          YM4 = YA(ICONA(4,IEL)) - YA(ICONA(1,IEL))
          ZM4 = ZA(ICONA(4,IEL)) - ZA(ICONA(1,IEL))
          V11 = XM1 - XM3
          V12 = YM1 - YM3
          V13 = ZM1 - ZM3
          V21 = XM2 - XM4
          V22 = YM2 - YM4
          V23 = ZM2 - ZM4
          SHLNRM(1) = V12*V23 - V22*V13
          SHLNRM(2) = V21*V13 - V11*V23
          SHLNRM(3) = V11*V22 - V21*V12
        END IF

C compute gradients for element IEL and put in SOLGRA

        CALL FLGRAD(IEL,ICOUNT,IELLST,CNTRA,SHLNRM,SOLEA,SOLGRA,
     &              ITT,iblk)

 110  CONTINUE
      RETURN
      END
