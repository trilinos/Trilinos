C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: vec.f,v $
C Revision 1.2  2009/03/25 12:36:49  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:17:27  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:59:29  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE VEC (IS3DIM, X0, Y0, Z0, XVAR, YVAR, ZVAR,
     &   VECSCL, VWSCL)
C=======================================================================

C   --*** VEC *** (DETOUR) Plot vector
C   --   Written by Amy Gilkey - revised 04/16/85
C   --
C   --VEC displays a vector starting at the given coordinates.  The vector
C   --represents the given X, Y, and Z values.
C   --
C   --Parameters:
C   --   IS3DIM - IN - true iff 3D versus 2D
C   --   X0, Y0, Z0 - IN - the vector coordinates
C   --   XVAR, YVAR, ZVAR - IN - the vector components
C   --   VECSCL - IN - the vector scale factor
C   --   VWSCL - IN - 1.0 if single view, 0.5 if multiple

      LOGICAL IS3DIM

      LOGICAL EXISTS

      EXISTS (M) = (MOD(M,2) .NE. 0)

      IF ((XVAR .EQ. 0.0) .AND. (YVAR .EQ. 0.0)) GOTO 100
      X1 = X0 + XVAR * VECSCL
      Y1 = Y0 + YVAR * VECSCL
      Z1 = Z0 + ZVAR * VECSCL
      CALL MP2VC (1, X0, Y0, X1, Y1, DX0, DY0, DX1, DY1, MASK)
      RAT = .0075 * VWSCL
      IF (EXISTS (MASK)) THEN
         CALL PLTARR (DX0, DY0, DX1, DY1, .5, RAT)
      END IF

  100 CONTINUE
      RETURN
      END
