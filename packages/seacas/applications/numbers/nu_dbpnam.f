C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBPNAM (OPTION, NVARGL, NVARNP, NVAREL,
     &                   NAMEGV, NAMENV, NAMEEV)
C=======================================================================

C   --*** DBPNAM *** (EXOLIB) Print database variable names
C   --   Written by Amy Gilkey - revised 01/21/88
C   --
C   --DBPNAM displays the database variable names.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options:
C   --                 'G' to print global variable names
C   --                 'N' to print nodal variable names
C   --                 'E' to print element variable names
C   --   NVARGL - IN - the number of global variables (if OPTION)
C   --   NVARNP - IN - the number of nodal variables (if OPTION)
C   --   NVAREL - IN - the number of element variables (if OPTION)
C   --   NAMEGV - IN - the global variable names (if OPTION)
C   --   NAMENV - IN - the nodal variable names (if OPTION)
C   --   NAMEEV - IN - the element variable names (if OPTION)

      include 'exodusII.inc'

      CHARACTER*(*) OPTION
      INTEGER NVARGL, NVARNP, NVAREL
      CHARACTER*(MXSTLN) NAMEGV(*)
      CHARACTER*(MXSTLN) NAMENV(*)
      CHARACTER*(MXSTLN) NAMEEV(*)
      LOGICAL ALL

      ALL = (OPTION .EQ. '*')

      WRITE (*, 10000)

      IF (ALL .OR. (INDEX (OPTION, 'G') .GT. 0)) THEN
         WRITE (*, 10010) 'Global: ', (NAMEGV(I), I=1,NVARGL)
      END IF
      IF (ALL .OR. (INDEX (OPTION, 'N') .GT. 0)) THEN
         WRITE (*, 10010) 'Nodal:  ', (NAMENV(I), I=1,NVARNP)
      END IF
      IF (ALL .OR. (INDEX (OPTION, 'E') .GT. 0)) THEN
         WRITE (*, 10010) 'Element:', (NAMEEV(I), I=1,NVAREL)
      END IF

      RETURN

10000  FORMAT (/, 1X, 'Variables Names:')

10010  FORMAT (4X, A, :, 2 (2X, A32), :, /,
     &        (12X, 2 (2X, A32)))
      END
