C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SHWINT (ITRANT, DIM3)
C=======================================================================

C   $Id: shwint.f,v 1.2 1991/01/09 12:59:38 gdsjaar Exp $
C   $Log: shwint.f,v $
C   Revision 1.2  1991/01/09 12:59:38  gdsjaar
C   Initial conversion from GEN3D to GENSHELL, no BC yet
C
c Revision 1.1.1.1  90/08/20  12:22:56  gdsjaar
c Gen3D Mesh Generation Program
c
c Revision 1.1  90/08/20  12:22:55  gdsjaar
c Initial revision
c

      CHARACTER*20 RSTR(9)
      CHARACTER*20 TYPE

      CALL NUMSTR (1, 4, DIM3, RSTR(1), LR)

      IF      (ITRANT .EQ.  0) THEN
         TYPE = 'Transform'
      ELSE IF (ITRANT .EQ.  1) THEN
         TYPE = 'Translate'
      ELSE IF (ITRANT .EQ.  2) THEN
         TYPE = 'Rotate'
      ELSE IF (ITRANT .EQ.  4) THEN
         TYPE = 'Warp'
      ELSE IF (ITRANT .EQ.  8) THEN
         TYPE = 'Twist'
      ELSE IF (ITRANT .EQ. 16) THEN
         TYPE = 'Project'
      ELSE IF (ITRANT .EQ. 32) THEN
         TYPE = 'ExpRotate'
      ELSE IF (ITRANT .EQ. 64) THEN
         TYPE = 'Spline'
      ELSE
         CALL PRTERR ('PROGRAM', 'Unknown transformation option')
         RETURN
      END IF

      LT = LENSTR(TYPE)

      WRITE (*, 20) TYPE(:LT), ' mesh, thickness = ', RSTR(1)(:LR)
 20   FORMAT (1X, 10A)
      RETURN
      END

