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

      SUBROUTINE BLDPTE(X,Y,Z,ICON,CENTER)
C
C***********************************************************************
C
C BLDPTE LOADS THE ELEMENT CENTROIDS INTO THE XYZSRF ARRAY
C
C Calls subroutine CNTR
C
C Called by MAPVAR
C
C***********************************************************************
C
      include 'bmesh.blk'
      include 'ebbyeb.blk'
      include 'tapes.blk'
C
      DIMENSION X(*),Y(*),Z(*),CENTER(NUMEBB,*),ICON(NELNDB,*)
      DIMENSION XX(27), YY(27), ZZ(27)
C
C
C      NNODES = NNELM(ITYPE)
      NNODES = NELNDB
      IF (ITYPE .EQ. 6) NNODES = 4
      DO 10 IEL = 1, NUMEBB
        DO 20 INOD = 1, NNODES
          XX(INOD) = X(ICON(INOD,IEL))
          YY(INOD) = Y(ICON(INOD,IEL))
          IF (NDIMB .EQ. 3)THEN
            ZZ(INOD) = Z(ICON(INOD,IEL))
          ELSE
            ZZ(INOD) = 0.
          END IF
   20   CONTINUE
        CALL CNTR(ITYPE,XX,YY,ZZ,CENTER(IEL,1),CENTER(IEL,2),
     &            CENTER(IEL,3))
   10 CONTINUE
      RETURN
      END
