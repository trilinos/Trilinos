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

      SUBROUTINE MKRNK(N,NTOTAL,NDIM,X,IND,IRNK,IRNK2)
C
C***********************************************************************
C
C     DESCRIPTION: THIS ROUTINE CONVERTS THE IRNK VECTORS IN THE SWEGLE
C                  SEARCH FROM AN INDIRECT GATHER TO A DIRECT GATHER
C
C       N        INTEGER   NUMBER OF ENTITIES THAT WAS SORTED
C       NDIM     INTEGER   NUMBER OF DIMENSIONS
C       X        REAL      ENTITIES TO BE SORTED
C
C       IND      INTEGER   INDEX VECTOR
C       IRNK     INTEGER   RANK VECTOR (INDIRECT)
C       IRNK2    INTEGER   RANK VECTOR (DIRECT)
C
C***********************************************************************
C
      include 'tapes.blk'
C
      DIMENSION X(NTOTAL,NDIM),IND(N,NDIM)
      DIMENSION IRNK(N,NDIM),IRNK2(N,NDIM,*)
C
       DO 11 IDM = 1, NDIM
         CALL INDEXX(N,X(1,IDM),IND(1,IDM),N)
         CALL RANK(N,IND(1,IDM),IRNK(1,IDM),N)
 11    CONTINUE
C
C CONSTRUCT DIRECT LISTS INTO ORDERED LIST OF POINTS
        IF(NDIM .EQ. 1)THEN
          DO 113 I = 1, N
            IRNK2(I,1,1) = IRNK(I,1)
113       CONTINUE
        ELSE IF( NDIM .EQ. 2)THEN
          DO 213 I = 1, N
            IRNK2(I,1,1) = IRNK(IND(I,1),2)
            IRNK2(I,2,1) = IRNK(IND(I,2),1)
213       CONTINUE
        ELSE IF( NDIM .EQ. 3)THEN
        DO 313 I=1,N
          IRNK2(I,1,1) = IRNK(IND(I,1),2)
          IRNK2(I,1,2) = IRNK(IND(I,1),3)
          IRNK2(I,2,1) = IRNK(IND(I,2),1)
          IRNK2(I,2,2) = IRNK(IND(I,2),3)
          IRNK2(I,3,1) = IRNK(IND(I,3),1)
          IRNK2(I,3,2) = IRNK(IND(I,3),2)
 313    CONTINUE
       ELSE
         PRINT*,'WRONG NUMBER OF DIMENSIONS IN MKRNK'
         STOP 
       ENDIF
C
      RETURN
      END
