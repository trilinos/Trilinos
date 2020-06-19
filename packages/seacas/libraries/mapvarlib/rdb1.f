C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

c
C======================================================================
*DECK,RDB1
      SUBROUTINE RDB1(XB,YB,ZB)
C
C     *****************************************************************
C
C     READS MESH-B
C
C     Called by MAPVAR
C
C     *****************************************************************
C
C     XB,etc  Original coordinates read from Mesh-B
C
C     *****************************************************************
C
      include 'ex2tp.blk'
C
      DIMENSION XB(*),YB(*),ZB(*)
C     *****************************************************************
C
C read coordinates
C
      CALL EXGCOR (NTP3EX,XB,YB,ZB,IERR)
C
      RETURN
      END
