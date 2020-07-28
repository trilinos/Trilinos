C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PRNAME (NOUT, NAMLEN,
     *                   NVARGL, NVARNP, NVAREL, NVARNS, NVARSS,
     &                   NAMEGV, NAMENV, NAMEEV, NAMNSV, NAMSSV)
C=======================================================================

C   --*** PRNAME *** (BLOT) Display database variable names
C   --   Written by Amy Gilkey - revised 01/14/88
C   --
C   --PRNAME displays the database variable names.
C   --
C   --Parameters:
C   --   NOUT - IN - the output file, <=0 for standard
C   --   NVARGL - IN - the number of global variables
C   --   NVARNP - IN - the number of nodal variables
C   --   NVAREL - IN - the number of element variables
C   --   NAMEGV - IN - the global variable names
C   --   NAMENV - IN - the nodal variable names
C   --   NAMEEV - IN - the element variable names

      CHARACTER*(NAMLEN) NAMEGV(*)
      CHARACTER*(NAMLEN) NAMENV(*)
      CHARACTER*(NAMLEN) NAMEEV(*)
      CHARACTER*(NAMLEN) NAMNSV(*)
      CHARACTER*(NAMLEN) NAMSSV(*)
      CHARACTER*128 FMT1, FMT

      IF (NOUT .GT. 0) WRITE (NOUT, 10000)

      IF (NOUT .GT. 0) THEN
         WRITE (NOUT, 10010)
      ELSE
         WRITE (*, 10010)
      END IF

      WRITE(FMT1,20) NAMLEN
      CALL SQZSTR(FMT1, LFMT)
      if (namlen .le. 20) then
         WRITE(FMT, 30) FMT1(:LFMT), FMT1(:LFMT)
      else
         WRITE(FMT, 40) FMT1(:LFMT), FMT1(:LFMT)
      endif

C ... Print them out.
      if (nout .le. 0) then
            WRITE (*, FMT) 'Global: ', (NAMEGV(I), I=1,NVARGL)
            WRITE (*, FMT) 'Nodal:  ', (NAMENV(I), I=1,NVARNP)
            WRITE (*, FMT) 'Element:', (NAMEEV(I), I=1,NVAREL)
            WRITE (*, FMT) 'Nodeset:', (NAMNSV(I), I=1,NVARNS)
            WRITE (*, FMT) 'Sideset:', (NAMSSV(I), I=1,NVARSS)
      else
            WRITE (NOUT, FMT) 'Global: ', (NAMEGV(I), I=1,NVARGL)
            WRITE (NOUT, FMT) 'Nodal:  ', (NAMENV(I), I=1,NVARNP)
            WRITE (NOUT, FMT) 'Element:', (NAMEEV(I), I=1,NVAREL)
            WRITE (NOUT, FMT) 'Nodeset:', (NAMNSV(I), I=1,NVARNS)
            WRITE (NOUT, FMT) 'Sideset:', (NAMSSV(I), I=1,NVARSS)
      end if

      RETURN

 20   FORMAT('A',I4)
 30   FORMAT ('(4X, A, :, 3 (2X, ',A,'), :, /,(12X, 3 (2X, ',A,')))')
 40   FORMAT ('(4X, A, :, 2 (2X, ',A,'), :, /,(12X, 2 (2X, ',A,')))')

10000  FORMAT (/, 1X, 'VARIABLES NAMES')
10010  FORMAT (/, 1X, 'Variables Names:')

      END
