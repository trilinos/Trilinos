C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      subroutine limits(STRING, ndim,
     &        xmin, xmax, ymin, ymax, zmin, zmax)
C=======================================================================
      character*(*) STRING

      WRITE (*, *) ' '
      WRITE (*, *) STRING(:LENSTR(STRING))
      WRITE (*, 10000) 'X', XMIN, 'X', XMAX, XMAX-XMIN
      WRITE (*, 10000) 'Y', YMIN, 'Y', YMAX, YMAX-YMIN
      IF (NDIM .EQ. 3) THEN
         WRITE (*, 10000) 'Z', ZMIN, 'Z', ZMAX, ZMAX-ZMIN
      END IF
10000 FORMAT( ' Minimum ',A1,' = ',1PE12.5,', Maximum ',A1,' = ',
     &     1PE12.5,', Range = ',1PE12.5)
      return
      end

