C    Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
C    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C    certain rights in this software
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C              
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C                            
C    * Neither the name of Sandia Corporation nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
C                                                    
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C    
      SUBROUTINE MYCOMP (MYV, VOID, LVOID,
     *   NVOIDS, DPOINT, LDICT, NNAMES, CHRCOL, LASTER)
C
      IMPLICIT INTEGER (A-Z)
      INCLUDE 'params.inc'
C
C     THIS ROUTINE PERFORMS THE NUMERIC DATA COMPRESSION OPERATION.
C
C************************************************************************
C
C     MYV      Reference array
C     VOID     Void table
C     LVOID    Dimension of VOID
C     NVOIDS   Number of voids
C     DPOINT   Dictionary pointer table
C     LDICT    Dimension of DPOINT
C     NNAMES   Number of names
C     CHRCOL   Column for character tables
C     LASTER   Error return code
C
      DIMENSION DPOINT(LDICT,CHRCOL,2), VOID(LVOID,CHRCOL,2)
      DIMENSION NNAMES(2), NVOIDS(2)
      CHARACTER*1 MYV(1)
C
C************************************************************************
C
      LASTER = SUCESS
C
C     The basic strategy is to look for an array in the dictionary
C     which is immediately preceeded by a void.  If found, a shift
C     is performed, and the void table is updated.
C
      IVOID = 0
  100 CONTINUE
      IVOID = IVOID + 1
  110 IF (IVOID .GT. NVOIDS(2)) GO TO 130
         VADDR = VOID(IVOID,2,1) + VOID(IVOID,2,2)
         DO 120 IDICT = 1, NNAMES(2)
            DADDR = DPOINT(IDICT,2,1)
            IF (VADDR .EQ. DADDR .AND. DPOINT(IDICT,2,2) .GT. 0) THEN
C
C              Perform data shift and update void table.
C
               CALL SHFTC (MYV, LDICT,
     *            DADDR, DADDR+DPOINT(IDICT,2,2)-1, VOID(IVOID,2,2))
               DPOINT(IDICT,2,1) = VOID(IVOID,2,1)
               VOID(IVOID,2,1) = DPOINT(IDICT,2,1) + DPOINT(IDICT,2,2)
               CALL VTABLE (0, 0, VOID(1,2,1), LVOID, NVOIDS(2),
     *            CHRCOL, LASTER)
               IF (LASTER .NE. SUCESS) RETURN
               GO TO 110
C
            END IF
  120    CONTINUE
      GO TO 100
  130 CONTINUE
      RETURN
      END
