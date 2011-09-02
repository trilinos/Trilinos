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
*DECK,WRTC
      SUBROUTINE WRTC (XB,YB,ZB,GVAR,SOLNB)
C     
C     ******************************************************************
C     
C     SUBROUTINE TO WRITE INTERPOLATED SOLUTION TO MESH-C EXODUS FILE
C
C     Called by MAPVAR
C    
C     ******************************************************************
C
C  GVAR   REAL Global variables
C  SOLNB  REAL Element variables interpolated onto Mesh-B
C  IXDIS  INT  Pointer to x-displ nodal variable
C  IYDIS  INT  Pointer to y-displ nodal variable
C  IZDIS  INT  Pointer to z-displ nodal variable
C
C     ******************************************************************
C
      include 'aexds1.blk'
      include 'bmesh.blk'
      include 'contrl.blk'
      include 'ex2tp.blk'
      include 'steps.blk'
      include 'varnpt.blk'
      include 'exodusII.inc'
C     
      DIMENSION SOLNB(NODESB,NVARNP)
      DIMENSION XB(*),YB(*),ZB(*),GVAR(*)
C     
C     ******************************************************************
      IF (ISTEP .EQ. -1)THEN
        NTM = NTIMES
      ELSE
        NTM = 1
      END IF
C
      DO 10 IST = 1, NTM
        IF (ISTEP .EQ. -1)THEN
          ISTP = IST
        ELSE
          ISTP = ISTEP
        END IF
c
c Time
c
        CALL EXGTIM (NTP2EX,ISTP,RTIME,IERR)
        IF (OUTTIM .LT. 0)THEN
          CALL EXPTIM (NTP4EX,IST,RTIME,IERR)
        ELSE
          CALL EXPTIM (NTP4EX,IST,OUTTIM,IERR)
        END IF
c
c Global variables
c         
        CALL EXGGV (NTP2EX,ISTP,NVARGP,GVAR,IERR)
        CALL EXPGV (NTP4EX,IST,NVARGP,GVAR,IERR)
C     
C Coordinates
C
        IF (IDEF .EQ. 1)THEN
          IF (IXDIS .NE. 0 .AND. IYDIS .NE. 0)THEN
            DO 60 I = 1, NODESB
              XB(I) = XB(I) - SOLNB(I,IXDIS)
              YB(I) = YB(I) - SOLNB(I,IYDIS)
              IF(NDIMB .EQ. 3) ZB(I) = ZB(I) - SOLNB(I,IZDIS)
 60         CONTINUE
          END IF
        END IF
        CALL EXPCOR(NTP4EX,XB,YB,ZB,IERR)
C     
 10   CONTINUE
C
       RETURN
       END
