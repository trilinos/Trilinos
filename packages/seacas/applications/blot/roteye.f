C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE ROTEYE (EYE, ROTCEN, ROTMAT, *)
C=======================================================================

C   --*** ROTEYE *** (MESH) Get rotation matrix from eye position
C   --   Written by Amy Gilkey - revised 10/21/86
C   --
C   --ROTEYE derives the rotation matrix from the given eye position.
C   --Rotation matrix gives an X rotation followed by a Y rotation.
C   --
C   --Parameters:
C   --   EYE - IN - the eye position
C   --   ROTCEN - IN - the center of rotation
C   --   ROTMAT - OUT - the rotation matrix
C   --   * - return statement if error in eye position

      REAL EYE(3)
      REAL ROTCEN(3)
      REAL ROTMAT(3,3)

      X = EYE(1) - ROTCEN(1)
      Y = EYE(2) - ROTCEN(2)
      Z = EYE(3) - ROTCEN(3)
      VMAG1 = SQRT (Y*Y + Z*Z)
      VMAG2 = SQRT (X*X + Y*Y + Z*Z)
      IF (VMAG1 .EQ. 0.0) GOTO 100
      COS1 = Z / VMAG1
      SIN1 = Y / VMAG1
      COS2 = VMAG1 / VMAG2
      SIN2 = -X / VMAG2

      ROTMAT(1,1) = COS2
      ROTMAT(2,1) = SIN1 * SIN2
      ROTMAT(3,1) = COS1 * SIN2
      ROTMAT(1,2) = 0.0
      ROTMAT(2,2) = COS1
      ROTMAT(3,2) = - SIN1
      ROTMAT(1,3) = - SIN2
      ROTMAT(2,3) = SIN1 * COS2
      ROTMAT(3,3) = COS1 * COS2

      RETURN

  100 CONTINUE
      CALL PRTERR ('CMDERR', 'Eye cannot be exactly on center')
      RETURN 1
      END
