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
*DECK, SETON1
      SUBROUTINE SETON1(CNTRA,SOLEA,SOLENA,IDBLK,XA,YA,ZA,ICONA,
     &                  NDLSTA,INVLN,INVCN,MAXLN,ISTP,ITT,iblk)
C
C  *********************************************************************
C
C  Subroutine SETON1 extracts nodal vaules of element variables by
C  performing a weighted least squares fit (4 or more elements) or 
C  a triangulation (3 elements) over the centroids of the elements 
C  attached to the current node.
C
C  Each element block must be processed independently in order to
C  avoid averaging element variables across material boundaries.
C  Note: the last set of DO loops acts over all nodes; to make sense
C        one element block must be completely processed before another
C        element block is sent into this subroutine.
C
C  Calls subroutines CNTR, VOL, EXTS, AVG, ERROR
C
C  Called by MAPVAR
C
C  *********************************************************************
C
C   CNTRA       a list of element centroid coordinates for all elements
C               in the current element block (1:ndima,1:numeba)
C   SOLEA       element variables (1:numeba,1:nvarel)
C   SOLENA      element variables at nodes (1:nodesa,1:nvarel)
C   IDBLK       current element block I.D.
C   XA,etc      coordinates
C   ICONA       connectivity array (1:nelnda,1:nodesa)
C   NDLSTA      list of nodes in element block - from RDA2 (1:numnda)
C   INVLN       number of elements connected to a node (1:numnda)
C   INVCN       inverse connectivity (1:maxln,1:numnda)
C   MAXLN       maximum number of elements connected to any node
C   ITT         truth table
C   iblk        element block being processed (not ID)
C
C** RELATIONSHIP BETWEEN NODAL IDENTIFICATIONS **
C   IGLND  = NDLSTA(INOD)   = ICONA(NOWLND,INVCN(1,IGLND))
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
      DIMENSION SOLENA(NODESA,NVAREL)
      DIMENSION XX(27), YY(27), ZZ(27)
      DIMENSION XA(*), YA(*), ZA(*), ICONA(NELNDA,*), NDLSTA(*)
      DIMENSION INVCN(MAXLN,*),INVLN(*), ITT(NVAREL,*)
C
C  *********************************************************************
C
      DO 10 I = 1, NODESA
      DO 10 J = 1, NVAREL
        SOLENA(I,J) = 0.
   10 CONTINUE
C
C  load up CNTRA array - coordinates of donor mesh element centroids
C
      NNODES = 4
        DO 20 IEL = 1, NUMEBA
          DO 25 I = 1, NNODES
            INODE = ICONA(I,IEL)
            XX(I) = XA(INODE)
            YY(I) = YA(INODE)
            ZZ(I) = ZA(INODE)
   25     CONTINUE
        CALL CNTR(13,XX,YY,ZZ,CNTRA(IEL,1),CNTRA(IEL,2),CNTRA(IEL,3))
   20   CONTINUE
C
C  put element variables into SOLEA array
C
      DO 30 IVAR = 1, NVAREL
        IF (ITT(IVAR,iblk) .EQ. 0)GO TO 30
        CALL EXGEV(NTP2EX,ISTP,IVAR,IDBLK,NUMEBA,SOLEA(1,IVAR),IERR)
C
        IF (NAMVAR(nvargp+IVAR)(1:6) .EQ. 'ELMASS') THEN
C
C replace element mass with density
C
          DO 35 IEL = 1, NUMEBA
            DO 40 I = 1, NNODES
              XX(I) = XA(ICONA(I,IEL))
              YY(I) = YA(ICONA(I,IEL))
              ZZ(I) = ZA(ICONA(I,IEL))
   40       CONTINUE
            CALL VOL(ITYPE,XX,YY,ZZ,VOLUME)
            SOLEA(IEL,IVAR) = SOLEA(IEL,IVAR) / VOLUME
   35     CONTINUE
        END IF
   30 CONTINUE
C
C start least squares extrapolation
C  Find the elements connected to the node. If fewer than 3 elements,
C  adjust search to find additional elements. If unable to get at 
C  least 3 elements, must be treated as special case (just average 
C  element values at node)(see below).
C
      DO 50 INOD = 1, NUMNDA
        IGLND = NDLSTA(INOD)
C
C  Process special case of only 1 element attached to node
C
        IF (INVLN(IGLND) .EQ. 1)THEN
C
C  Get node number diagonally across element, in most cases this
C  node will have 4 elements attached.
C
          DO 70 I = 1, NNODES
            IF (IGLND .EQ. ICONA(I,INVCN(1,IGLND)))NXTLND = I + 2
   70     CONTINUE
          IF (NXTLND .GT. NNODES) NXTLND = NXTLND - NNODES
          NXGLND = ICONA(NXTLND,INVCN(1,IGLND))
C
C  If 3 or more elements,
c  fit a plane through the element centroids, project element
c  centroids and original node onto plane, extrapolate 
c  in 2-d (coords of plane) to original node (done in EXTS).
c  If 2 or less elements, 
c  average original element variables at original node
c
          IF (INVLN(NXGLND) .GT. 2)THEN
            CALL EXTS(IGLND,INVCN,MAXLN,NXGLND,INVLN(NXGLND),XA,YA,ZA,
     &                CNTRA,SOLEA,SOLENA,ITT,iblk)
          ELSE
            CALL AVG(IGLND,INVCN,MAXLN,INVLN(IGLND),SOLEA,SOLENA,
     &               ITT,iblk)
          END IF
C
C  Process special case of only 2 elements attached to node
C
        ELSE IF (INVLN(IGLND) .EQ. 2)THEN
c
c  get second node that is shared by both elements. That is the
c  node on the other end of the shared element side.
c
          DO 90 I = 1, NNODES
          DO 90 J = 1, NNODES
            IF(ICONA(I,INVCN(1,IGLND)) .NE. IGLND .AND. 
     &         ICONA(I,INVCN(1,IGLND)) .EQ. ICONA(J,INVCN(2,IGLND)))
     &         NXGLND = ICONA(I,INVCN(1,IGLND))
   90     CONTINUE
c
c  If this second node has more than 3 elements, extrapolate. Otherwise
c  average. (at original node)
c
          IF (INVLN(NXGLND) .GT. 2)THEN
            CALL EXTS(IGLND,INVCN,MAXLN,NXGLND,INVLN(NXGLND),XA,YA,ZA,
     &                CNTRA,SOLEA,SOLENA,ITT,iblk)
          ELSE
            CALL AVG(IGLND,INVCN,MAXLN,INVLN(IGLND),SOLEA,SOLENA,
     &               ITT,iblk)
          END IF
        ELSE
          CALL EXTS(IGLND,INVCN,MAXLN,IGLND,INVLN(IGLND),XA,YA,ZA,
     &              CNTRA,SOLEA,SOLENA,ITT,iblk)
        END IF
   50 CONTINUE
      RETURN
      END
