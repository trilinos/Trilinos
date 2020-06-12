C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Log: grbox.f,v $
C Revision 1.2  2009/03/25 12:36:44  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.1  1994/04/07 20:02:06  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:51:26  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE GRBOX (LINSOL, XMIN, XMAX, YMIN, YMAX)
C=======================================================================

C   --*** GRBOX *** (GRPLIB) Draw box (either outline or solid)
C   --   Written by Amy Gilkey - revised 02/20/87
C   --
C   --GRBOX draws a box (either the outline or a solid) defined by two
C   --coordinate pairs.
C   --
C   --Parameters:
C   --   LINSOL - IN - 'L' if outline, 'S' if solid
C   --   XMIN, XMAX - IN - the horizontal minimum and maximum of box
C   --      (in device coordinates)
C   --   YMIN, YMAX - IN - the vertical minimum and maximum of box
C   --      (in device coordinates)

C   --Routines Called:
C   --   PLTPLY - (PLTLIB) Plot a filled polygon
C   --   PLTVCT - (PLTLIB) Plot a vector

      CHARACTER LINSOL
      REAL XMIN, XMAX, YMIN, YMAX

      REAL X(5), Y(5)

      X(1) = XMIN
      Y(1) = YMIN
      X(2) = XMAX
      Y(2) = Y(1)
      X(3) = X(2)
      Y(3) = YMAX
      X(4) = X(1)
      Y(4) = Y(3)
      IF (LINSOL .NE. 'S') THEN
         X(5) = X(1)
         Y(5) = Y(1)
      END IF

      IF (LINSOL .NE. 'S') THEN
         CALL PLTVCT (4, X(1), Y(1), X(2), Y(2))
      ELSE
         CALL PLTPLY (4, X(1), Y(1), X(2), Y(2))
      END IF

      RETURN
      END
