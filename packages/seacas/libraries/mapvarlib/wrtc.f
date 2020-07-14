C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

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
        if (nvargp .gt. 0) then
           CALL EXGGV (NTP2EX,ISTP,NVARGP,GVAR,IERR)
           CALL EXPGV (NTP4EX,IST,NVARGP,GVAR,IERR)
        end if
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
