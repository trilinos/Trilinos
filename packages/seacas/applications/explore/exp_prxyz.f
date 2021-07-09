C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE PRXYZ (OPTION, NOUT, NDIM, NAMECO, NUMNP, LISNP, CORD,
     *  MAP, DOMAP)
C=======================================================================

C   --*** PRXYZ *** (EXPLORE) Display database coordinates
C   --
C   --PRXYZ displays the coordinate array.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options
C   --   NOUT - IN - the output file, <=0 for standard
C   --   NDIM - IN - the number of coordinates per node
C   --   NAMECO - IN - the coordinate names
C   --   NUMNP - IN - the number of nodes
C   --   LISNP - IN - the indices of the selected nodes
C   --   CORD - IN - the nodal coordinates

      include 'exodusII.inc'
      CHARACTER*(*) OPTION
      CHARACTER*(*) NAMECO(*)
      INTEGER LISNP(0:*)
      REAL CORD(NUMNP,NDIM)
      INTEGER MAP(*)
      LOGICAL DOMAP
      INTEGER GETPRC, PRTLEN
      CHARACTER*128 FMT1, FMT

      PRTLEN = GETPRC() + 7
      WRITE(FMT1,20) PRTLEN, PRTLEN-7
      CALL SQZSTR(FMT1, LFMT)
      WRITE(FMT, 30) FMT1(:LFMT), FMT1(:LFMT)

      IF (NOUT .GT. 0) WRITE (NOUT, 10000)

      if (domap) then
        if (nout .gt. 0) then
          write (nout, 10005)
        else
          write (*, 10005)
        end if
      end if
      IF (NOUT .GT. 0) THEN
         WRITE (NOUT, 10010) (NAMECO(I)(:8), I=1,NDIM)
      ELSE
         WRITE (*, 10010) (NAMECO(I)(:8), I=1,NDIM)
      END IF

      DO 100 IX = 1, LISNP(0)
         INP = LISNP(IX)
         if (domap) then
           id = map(inp)
         else
           id = inp
         end if
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, FMT, IOSTAT=IDUM)
     &         ID, (CORD(INP,I), I=1,NDIM)
         ELSE
            WRITE (*, FMT, IOSTAT=IDUM)
     &         ID, (CORD(INP,I), I=1,NDIM)
         END IF
  100 CONTINUE

      RETURN

10000  FORMAT (/, 1X, 'COORDINATES')
10005  FORMAT (1X, 'Nodal ids are Global')
 20   FORMAT('1PE',I2.2,'.',I2.2)
 30   FORMAT('(1X, ''Node'', I12, 5 (2X, ',A,'), :, /,',
     $     '(15X, 5 (2X, ',A,')))')

10010  FORMAT (/, 1X, 4X, 5X, 4X, 5 (2X, A8, :, 7X), :, /,
     &   (1X, 6 (7X, A8, :, 9X)))
      END
