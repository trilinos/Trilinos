C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
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
C 

C=======================================================================
*DECK, ELGRAD
      SUBROUTINE ELGRAD(CNTRA,XA,YA,ZA,SOLEA,SOLGRA,ICHKEL,IDBLK,
     &                  ICONA,INVLN,INVCN,MAXLN,ISTP,ITT,IM)
C
C  *********************************************************************
C
C  Subroutine ELGRAD computes the coefficients of the gradient for
C  each variable within an element based on a constrained (to the
C  value of that element at its centroid) least squares fit to the
C  values in the surrounding elements. If not enough data is available
C  to perform the least squares is a direction, that coefficient is set
C  to zero.
C
C  Calls subroutines CNTR, VOL, FLGRAD, ERROR
C
C  Called by MAPVAR
C
C  *********************************************************************
C
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
C   IM        element block being processed (not ID)
C
C
C  *********************************************************************
C
      include 'exodusII.inc'
C
      include 'aexds1.blk'
      include 'aexds2.blk'
      include 'amesh.blk'
      include 'ebbyeb.blk'
      include 'ex2tp.blk'
      include 'tapes.blk'
C
      DIMENSION CNTRA(NUMEBA,*), SOLEA(NUMEBA,*)
      DIMENSION SOLGRA(NDIMA,NUMEBA,*), ICHKEL(*)
      DIMENSION XX(27), YY(27), ZZ(27), IELLST(100), SHLNRM(3)
      DIMENSION XA(*), YA(*), ZA(*), ICONA(NELNDA,*)
      DIMENSION INVCN(MAXLN,*),INVLN(*), ITT(NVAREL,*)
C
C  *********************************************************************
C
      IF (ITYPE .EQ. 4 .OR. ITYPE .EQ. 5)THEN
        CALL ERROR('ELGRAD','ELEMENT TYPE',' ',ITYPE,
     &             'ELEMENT VARIABLE PROCESSING NOT YET IMPLEMENTED',
     &              0,' ',' ',1)
      END IF
C
C initialize ICHKEL array - only needs to be done once because 
C checkmark used is specfic to each element
C
      DO 10 I = 1, NUMEBA
        ICHKEL(I) = 0
 10   CONTINUE
C
C  load up CNTRA array - coordinates of mesh-A element centroids
C
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
C
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
C
   40   CONTINUE
      END IF
C
C  put element variables into SOLEA array
C
      DO 80 IVAR = 1, NVAREL
        IF (ITT(IVAR,IM) .EQ. 0)GO TO 80
        CALL EXGEV(NTP2EX,ISTP,IVAR,IDBLK,NUMEBA,SOLEA(1,IVAR),IERR)
C
        IF (NAMVAR(nvargp+IVAR)(1:6) .EQ. 'ELMASS') THEN
C
C replace element mass with density
C
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
C
            SOLEA(IEL,IVAR) = SOLEA(IEL,IVAR) / VOLUME
   90     CONTINUE
        END IF
   80 CONTINUE
C
C start element centroid based gradient calculation 
C 
C get list of elements that share at least one node with current
C element being processed IELLST
C
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
C
C If shells element, compute element normal
C use cross product of vectors from mid-sides
C
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
C
C compute gradients for element IEL and put in SOLGRA
C
        CALL FLGRAD(IEL,ICOUNT,IELLST,CNTRA,SHLNRM,SOLEA,SOLGRA,
     &              ITT,IM)
C
 110  CONTINUE
      RETURN
      END
