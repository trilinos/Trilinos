C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE PRNODE (OPTION, NOUT, NUMNP, LISNP,
     &   NVARNP, LISNV, NAMENV, VARNP, MAP, DOMAP)
C=======================================================================

C   --*** PRNODE *** (EXPLORE) Display current database nodal variables
C   --
C   --PRNODE displays the nodal data for a time step.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options
C   --   NOUT - IN - the output file, <=0 for standard
C   --   NUMNP - IN - the number of nodes
C   --   LISNP - IN - the indices of the selected nodes
C   --   NVARNP - IN - the number of nodal variables
C   --   LISNV - IN - the indices of the selected nodal variables
C   --   NAMENV - IN - the names of the nodal variables
C   --   VARNP - IN - the nodal variables for the time step

      include 'exodusII.inc'
      CHARACTER*(*) OPTION
      INTEGER LISNP(0:*)
      INTEGER LISNV(0:*)
      CHARACTER*(*) NAMENV(*)
      REAL VARNP(numnp,*)
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

      do 90 i=1, lisnv(0)
        irow = ((i-1)/5)+1
        icol = i - (irow-1)*5
        IF (NOUT .GT. 0) THEN
           WRITE (NOUT, 10010) irow, icol, NAMENV(LISNV(I))
        ELSE
           WRITE (*, 10010) irow, icol, NAMENV(LISNV(I))
        END IF
 90   continue

      DO 100 IX = 1, LISNP(0)
         INP = LISNP(IX)
         if (domap) then
            id = map(inp)
         else
            id = inp
         end if

         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, FMT, IOSTAT=IDUM)
     &           ID, (VARNP(INP, LISNV(I)), I=1,LISNV(0))
         ELSE
            WRITE (*, FMT, IOSTAT=IDUM)
     &           ID, (VARNP(INP, LISNV(I)), I=1,LISNV(0))
         END IF
 100  CONTINUE

      RETURN

 20   FORMAT('1PE',I2.2,'.',I2.2)
 30   FORMAT('(1X, ''Node'', I12, 5 (2X, ',A,'), :, /,',
     $     '(15X, 5 (2X, ',A,')))')
10000 FORMAT (/, 1X, 'NODAL TIME STEP VARIABLES')
10005 FORMAT (1X, 'Nodal ids are Global')
10010 FORMAT (1X, 'Row ',I4,', Column ',I1,' is variable ',A)
      END
