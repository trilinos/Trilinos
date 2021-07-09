C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE PRNAME (OPTION, NOUT,
     &     NAMEGV, NAMENV, NAMEEV, NAMENS, NAMESS)
C=======================================================================

C   --*** PRNAME *** (EXPLORE) Display database variable names
C   --
C   --PRNAME displays the database variable names.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options:
C   --      'G' to print global variable names
C   --      'N' to print nodal variable names
C   --      'E' to print element variable names
C   --   NOUT - IN - the output file, <=0 for standard
C   --   NVARGL - IN - the number of global variables
C   --   NVARNP - IN - the number of nodal variables
C   --   NVAREL - IN - the number of element variables
C   --   NAMEGV - IN - the global variable names
C   --   NAMENV - IN - the nodal variable names
C   --   NAMEEV - IN - the element variable names

      include 'exp_dbnums.blk'

      CHARACTER*(*) OPTION
      CHARACTER*(*) NAMEGV(*)
      CHARACTER*(*) NAMENV(*)
      CHARACTER*(*) NAMEEV(*)
      CHARACTER*(*) NAMENS(*)
      CHARACTER*(*) NAMESS(*)
      CHARACTER*128 FMT1, FMT

      IF (NOUT .GT. 0) WRITE (NOUT, 10000)

      IF (NOUT .GT. 0) THEN
         WRITE (NOUT, 10010)
      ELSE
         WRITE (*, 10010)
      END IF

      WRITE(FMT1,20) NAMLEN
      CALL SQZSTR(FMT1, LFMT)
      WRITE(FMT, 30) FMT1(:LFMT), FMT1(:LFMT)

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'G') .GT. 0)) THEN
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, FMT) 'Global: ', (NAMEGV(I), I=1,NVARGL)
         ELSE
            WRITE (*, FMT) 'Global: ', (NAMEGV(I), I=1,NVARGL)
         END IF
      END IF
      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'N') .GT. 0)) THEN
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, FMT) 'Nodal:  ', (NAMENV(I), I=1,NVARNP)
         ELSE
            WRITE (*, FMT) 'Nodal:  ', (NAMENV(I), I=1,NVARNP)
         END IF
      END IF
      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'E') .GT. 0)) THEN
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, FMT) 'Element:', (NAMEEV(I), I=1,NVAREL)
         ELSE
            WRITE (*, FMT) 'Element:', (NAMEEV(I), I=1,NVAREL)
         END IF
      END IF
      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'M') .GT. 0)) THEN
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, FMT) 'Nodeset:', (NAMENS(I), I=1,NVARNS)
         ELSE
            WRITE (*, FMT) 'Nodeset:', (NAMENS(I), I=1,NVARNS)
         END IF
      END IF
      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'S') .GT. 0)) THEN
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, FMT) 'Sideset:', (NAMESS(I), I=1,NVARSS)
         ELSE
            WRITE (*, FMT) 'Sideset:', (NAMESS(I), I=1,NVARSS)
         END IF
      END IF

      RETURN

 20   FORMAT('A',I4)
 30   FORMAT ('(4X, A, :, 2 (2X, ',A,'), :, /,(12X, 2 (2X, ',A,')))')

10000  FORMAT (/, 1X, 'VARIABLES NAMES')
10010  FORMAT (/, 1X, 'Variables Names:')
      END
