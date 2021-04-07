C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE PRFRM (NOUT)
C=======================================================================

C   --*** PRFRM *** (EXPLORE) Display Coordinate Frame information
C   --
C   --Parameters:
C   --   NOUT - IN - the output file, <=0 for standard

      include 'exp_dbase.blk'
      include 'exodusII.inc'

      REAL RDUM
      character*1 cdum

C ... Get count of coordinate frames in model...
      call exinq(ndb, EXNCF, ncf, rdum, cdum, ierr)

      call prfrm1(nout, ncf)
      end

C=======================================================================
      subroutine prfrm1(nout, ncf)
C=======================================================================
      include 'exp_dbase.blk'
      include 'exodusII.inc'

      integer cfids(ncf), tags(ncf)
      real    coord(27)

      INTEGER GETPRC, PRTLEN
      CHARACTER*256 FMT1, FMT

      character*12 tag
      character*32 str32

      PRTLEN = GETPRC() + 7
      WRITE(FMT1,20) PRTLEN, PRTLEN-7
      CALL SQZSTR(FMT1, LFMT)
      WRITE(FMT, 35) 'I12', FMT1(:LFMT), FMT1(:LFMT), FMT1(:LFMT)

      if (nout .gt. 0) then
        WRITE (nout, 10000)
      end if

      call exgfrm(ndb, ncf, cfids, coord, tags, ierr)

      CALL INTSTR (1, 0, NCF, STR32, LSTR)
      IF (NOUT .GT. 0) THEN
        WRITE (NOUT, 10010) STR32(:LSTR)
      ELSE
        WRITE (*, 10010) STR32(:LSTR)
      END IF

      do i=1, ncf
        icbeg = 9*(i-1)+1
        icend = 9*i
        if (tags(i) .eq. EXCFREC) then
          tag = 'Rectangular'
        else if (tags(i) .eq. EXCFCYL) then
          tag = 'Cylindrical'
        else if (tags(i) .eq. EXCFSPH) then
          tag = 'Spherical  '
        end if

        IF (NOUT .GT. 0) THEN
          write (nout,FMT) cfids(i), tag, (coord(j),j=icbeg, icend)
        ELSE
          write (*,FMT) cfids(i), tag, (coord(j),j=icbeg, icend)
        END IF
      end do

      RETURN

 20   FORMAT('1PE',I2.2,'.',I2.2)

 35   FORMAT ('(/,'' Coordinate Frame '',',A,
     *  ','': '',A,/,5x,''Origin:          '',3(1x, ',A,
     *  '),/,5x,''3rd Axis Point:  '',3(1x, ',A,
     *  '),/,5x,''1-3 Plane Point: '',3(1x, ',A,'))')

10000  FORMAT (/, 1X, 'COORDINATE FRAMES')
10010  FORMAT (' Number of Coordinate Frames = ', A)
      END
