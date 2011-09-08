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

C========================================================================
*DECK, ELTON0
      SUBROUTINE ELTON0(ICONA,NELTN,SOLEA,SOLENA,
     &                  IDBLK,XA,YA,ZA,ISTP,ITT,IM)
C
C  *********************************************************************
C
C  Subroutine ELTON0 extracts nodal vaules of element variables by
C  looping over each element and summing the value of the variable
C  at that element to each node in the connectivity list for that 
C  element. Then the nodal summation of element variables is divided
C  by the number of elements that contributed to that node (resulting
C  in a nodal average of the element value.) This is done for the old
C  mesh elements and nodes to facilitate interpolation.
C
C  Each element block must be processed independently in order to
C  avoid averaging element variables across material boundaries.
C  Note: the last set of DO loops acts over all nodes; to make sense
C        one element block must be completely processed before another
C        element block is sent into this subroutine.
C
C  Calls subroutine VOL
C
c  Called by MAPVAR
C
C  *********************************************************************
C
C   ICONA       mesh-A connectivity (1:nelnda,1:numeba)
C   NELTN       number of elements tied to each node (1:nodesa)
C   SOLEA       element variables (1:numeba,1:nvarel)
C   SOLENA      element variables at nodes (1:nodesa,1:nvarel)
C   IDBLK       current element block I.D.
C   XA,YA,ZA    coordinates
C   XX,YY,ZZ    vector of coordinates of nodes for an element
C   ISTP        current time step
C   ITT         truth table
C   IM          element block being processed (not ID)
C
C  *********************************************************************
C
      include 'exodusII.inc'

      include 'aexds1.blk'
      include 'aexds2.blk'
      include 'amesh.blk'
      include 'ebbyeb.blk'
      include 'ex2tp.blk'
      include 'tapes.blk'
C
      DIMENSION ICONA(NELNDA,*), NELTN(*)
      DIMENSION SOLEA(NUMEBA,*), SOLENA(NODESA,NVAREL), ITT(NVAREL,*)
      DIMENSION XA(*), YA(*), ZA(*), XX(27), YY(27), ZZ(27)
C
C  *********************************************************************
C
      IF (ITYPE .EQ. 4 .OR. ITYPE .EQ. 5)THEN
        CALL ERROR('ELTON0','ELEMENT TYPE',' ',ITYPE,
     &             'ELEMENT VARIABLE PROCESSING NOT YET IMPLEMENTED',
     &              0,' ',' ',1)
      END IF
C
C
      DO 10 I = 1, NODESA
        NELTN(I) = 0
      DO 10 J = 1, NVAREL
        SOLENA(I,J) = 0.
   10 CONTINUE
C
C      NNODES = NNELM(ITYPE)
      NNODES = NELNDA
      IF (ITYPE .EQ. 6) NNODES = 4
      DO 20 NEL = 1, NUMEBA
        DO 30 I = 1, NNODES
C
C  number of elements associated with each node - used for
C    computing an average later on
C
          NELTN(ICONA(I,NEL)) = NELTN(ICONA(I,NEL)) + 1
   30   CONTINUE
   20 CONTINUE
C
      DO 40 IVAR = 1, NVAREL
        IF (ITT(IVAR,IM) .EQ. 0)GO TO 40
        CALL EXGEV(NTP2EX,ISTP,IVAR,IDBLK,NUMEBA,SOLEA(1,IVAR),IERR)
C
        IF (NAMVAR(nvargp+IVAR)(1:6) .EQ. 'ELMASS') THEN
C
C replace element mass with nodal density for interpolation
C
          DO 50 IEL = 1, NUMEBA
            DO 60 I = 1, NNODES
              XX(I) = XA(ICONA(I,IEL))
              YY(I) = YA(ICONA(I,IEL))
              IF (NDIMA .EQ. 3)THEN
                ZZ(I) = ZA(ICONA(I,IEL))
              ELSE
                ZZ(I) = 0.
              END IF
   60       CONTINUE
            CALL VOL(ITYPE,XX,YY,ZZ,VOLUME)
            SOLEA(IEL,IVAR) = SOLEA(IEL,IVAR) / VOLUME
   50     CONTINUE
        END IF
C
C  accumulate element variables to nodes
C
        DO 70 NEL = 1, NUMEBA
          DO 80 I = 1, NNODES
            SOLENA(ICONA(I,NEL),IVAR) = 
     &      SOLENA(ICONA(I,NEL),IVAR) + SOLEA(NEL,IVAR)
   80     CONTINUE
   70   CONTINUE
C
C  divide by number of elements contributing to each node (average)
C
        DO 90 I = 1, NODESA
          IF(NELTN(I) .NE. 0)THEN
            SOLENA(I,IVAR) = SOLENA(I,IVAR) / float(NELTN(I))
          END IF
   90   CONTINUE
   40 CONTINUE
      RETURN
      END

