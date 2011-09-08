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
      SUBROUTINE MXGET (MYLOC, MNGET, VOID, LVOID, NVOIDS,
     *   CHRCOL, LASTER, VROW)
C
      IMPLICIT INTEGER (A-Z)
      INCLUDE 'params.inc'
C
C     This subroutine returns the location (row number) of a void with
C     sufficient space for the memory request.  If necessary, memory is
C     acquired from the system.  The memory is contiguous.
C
C
C***********************************************************************
C
C     MYLOC    Address of internal reference array
C     MNGET    Memory request in numerical storage units
C     VOID     Void table
C     LVOID    Dimension of void table
C     NVOIDS   Number of voids
               DIMENSION VOID(LVOID,CHRCOL,2), NVOIDS(2)
C     CHRCOL   Column for character tables.
C     LASTER   Error return
C     VROW     Row number of void which satisfies the memory request
C
C***********************************************************************
C
C     IS THE MEMORY REQUEST SENSIBLE?
C
      IF (MNGET .LT. 0) THEN
         LASTER = BADLEN
         RETURN
      ELSE IF (MNGET .EQ. 0) THEN
         LASTER = SUCESS
         RETURN
      END IF
C
      CALL MXLOOK (MNGET, VOID, CHRCOL*LVOID, NVOIDS(1), VROW, LASTER)
      IF (LASTER .EQ. SUCESS) RETURN
C
C     CALL EXTENSION LIBRARY ROUTINE TO GET SPACE FROM SYSTEM.
C
      CALL EXMEMY (MNGET, LOC, MEMRET)
      LOC = LOC - MYLOC + 1
C
c  On return from exmemy, memret is set equal to -1 on an invalid
c  memory request (at least that's the plan under the new C code
c  extension library).  Therefore, I've made the change that should
c  test the appropriate condition.
c
      IF (MEMRET .LT. 0) THEN
C
C        ILLEGAL MEMORY REQUEST.
C
         LASTER = NOGET
         RETURN
C
      END IF
C
C     UPDATE VOID TABLE.
C
      CALL VTABLE (LOC, MEMRET, VOID, LVOID, NVOIDS(1), CHRCOL, LASTER)
      IF (LASTER .NE. SUCESS) RETURN
C
      CALL MXLOOK (MNGET, VOID, CHRCOL*LVOID, NVOIDS(1), VROW, LASTER)
C
      RETURN
      END
