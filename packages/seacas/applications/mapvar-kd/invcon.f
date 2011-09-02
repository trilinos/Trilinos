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
      SUBROUTINE INVCON(INVLN,MAXLN,INVCN,ICONA)
C
C************************************************************************
C
C Subroutine INVC0N computes the inverse connectivity (elements connected
C to a node). 
C
c Called by MAPVAR
C
C Calls ERROR
C
C************************************************************************
C
C  INVLN  INT   The number of elements connected to a node (1:numnda)
C  MAXLN  INT   The maximum number of elements connected to any node
C  INVCN  INT   The inverse connectivity (1:maxln,1:numnda)
C  ICONA  INT   The connectivity array (1:nelnda,1:numela)
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
      DO 5 I = 1, NODESA
        INVLN(I) = 0
      DO 5 J = 1, MAXLN
        INVCN(J,I) = 0
    5 CONTINUE
C
C      NNODES = NNELM(ITYPE)
      NNODES = NELNDA
      IF (ITYPE .EQ. 6) NNODES = 4
      DO 20 J = 1, NUMEBA
        DO 10 I = 1, NNODES
          node = icona(i,j) 
          IF (invln(node) .eq. 0 .or.
     *      INVCN(INVLN(node),node) .NE. J) THEN
            INVLN(node) = INVLN(node) + 1
            IF (INVLN(node) .GT. MAXLN) 
     &        CALL ERROR('INVCON',' ',
     &        'TOO MANY ELEMENTS CONNECTED TO NODE',
     &        node,'INVCN ARRAY DIMENSIONED FOR NO MORE THAN',
     &        MAXLN,'RESET IN SUBROUTINE RDA2',' ',1)
            INVCN(INVLN(node),node) = J
          END IF
 10     CONTINUE
 20   CONTINUE
      RETURN
      END
      
