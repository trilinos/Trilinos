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
      SUBROUTINE AVG(IGLND,INVCN,MAXLN,INVLEN,SOLEA,SOLENA,ITT,iblk)
C
C************************************************************************
C
C Subroutine AVG provides for translating nodal values of element
C variables back to the element centroids for the special case where 
C too few elements can be associated with a node. Element variable 
C data is simply averaged at that node.
C
C Called by ELTON1
C
C************************************************************************
C
C  IGLND  INT   The global node number
C  INVCN  INT   The inverse connectivity (1:maxln,1:numnda)
C  MAXLN  INT   The maximum nomber of elements connected to any node
C  INVLEN INT   The number of elements connected to this node
C  SOLEA  REAL  Element variables (1:numeba,1:nvarel)
C  SOLENA REAL  Element variables at nodes" (1:nodesa,1:nvarel)
C  NDLSTA INT   The array that identifies the local element block node
C               number with the global mesh node number (1:numnda)
C  ITT    INT   Truth table
C  iblk   INT   Block number being processed (not block ID)
C
C************************************************************************
C
      include 'aexds1.blk'
      include 'amesh.blk'
      include 'ebbyeb.blk'
C
      DIMENSION INVCN(MAXLN,*),SOLEA(NUMEBA,*),
     &          SOLENA(NODESA,NVAREL), ITT(NVAREL,*)
C
C************************************************************************
C
      DO 10 IVAR = 1, NVAREL
        IF (ITT(IVAR,iblk) .EQ. 0)GO TO 10
        SUM = 0.
        DO 20 J = 1, INVLEN
          SUM = SUM + SOLEA(INVCN(J,IGLND),IVAR)
   20   CONTINUE
        SOLENA(IGLND,IVAR) = SUM / INVLEN
   10 CONTINUE
      RETURN
      END
