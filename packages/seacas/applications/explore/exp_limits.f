C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      subroutine limits(ndim, numnp, cord)
C=======================================================================
      REAL CORD(NUMNP, NDIM)

      call minmax(numnp, cord(1,1), xmin, xmax)
      call minmax(numnp, cord(1,2), ymin, ymax)
      if (ndim .eq. 3) then
        call minmax(numnp, cord(1,3), zmin, zmax)
      else
        zmin = 0
        zmax = 0
      end if

      WRITE (*, *) ' '
      WRITE (*, *) 'Mesh Limits: '
      WRITE (*, 10000) 'X', XMIN, 'X', XMAX, XMAX-XMIN
      IF (NDIM .GE. 2) THEN
        WRITE (*, 10000) 'Y', YMIN, 'Y', YMAX, YMAX-YMIN
      ENDIF
      IF (NDIM .EQ. 3) THEN
         WRITE (*, 10000) 'Z', ZMIN, 'Z', ZMAX, ZMAX-ZMIN
      END IF
10000 FORMAT( ' Min ',A1,' = ',1PE16.9,', Max ',A1,' = ',
     &     1PE16.9,', Range = ',1PE16.9)
      return
      end

