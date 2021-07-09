C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE WRTCSV (NCRV, MAXPTS, NPTS, PTS, TXLAB, NAMES, DOLEGN,
     *  MAPEL, MAPND)
C=======================================================================

C   --*** WRTCSV *** (XYPLOT) Write curve to CSV neutral file
C   --
C   --WRTCSV writes the data for a curve to a neutral file which is
C   --in a comma-separated format. The first
C   --time the routine is called, the neutral file is opened.
C   --
C   --Parameters:
C   --   NPTS - IN - the number of points on the curve
C   --   PTS - IN - the plot data;
C   --      PTS(x,NTPVAR+1) holds the times if TIMPLT
C   --      PTS(x,NTPVAR+2) holds the compressed times if TIMPLT and needed
C   --   PLTITL - IN - the plot title describing the curve
C   --      (e.g. "TIME vs SIGXX at ELEMENT 30")
C   --   TXLAB, TYLAB - IN - the X and Y axis labels, either the
C   --      user-input labels or the plot variable descriptions
C   --

      PARAMETER (MAXLEG = 4096)
      include 'dbname.blk'
      include 'dbtitl.blk'
      include 'legopt.blk'
      include 'xylim.blk'
      include 'csv.blk'

      include 'params.blk'
      include 'tpvars.blk'

      CHARACTER*2048 filnam, errmsg
      INTEGER NPTS(*)
      REAL PTS(MAXPTS,*)
      CHARACTER*(*) TXLAB, NAMES(*)
      LOGICAL DOLEGN
      INTEGER MAPEL(*), MAPND(*)

      CHARACTER*(MAXLEG) LEGEND
      CHARACTER*(MXLNLN) PV

      IF (.NOT. CSVOPN) THEN
        filnam = basenam(:lenstr(basenam)) // '.csv'

C      --Open the csv neutral file and write the title line

        write (*,*) "CSV File: ", filnam(:lenstr(filnam))
        open (unit=ncsv, file=filnam(:lenstr(filnam)), form='formatted',
     *    status='unknown', iostat=ierr)
        IF (IERR .NE. 0) THEN
          ERRMSG = 'Neutral CSV file "'//FILNAM(:LENSTR(FILNAM))//
     *      '" could not be opened.'
          CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
          GOTO 170
        END IF
        CSVOPN = .TRUE.

      END IF

      legend = 'TIME'
      if (timplt) then
        n = 1
      else
        n = 2
      end if

      do i = n, ncrv+n-1
        call tplabv(0, itvid(i), names(itvid(i)), itvne(i), pv,
     *    MAPEL, MAPND)
        lpv = lenstr(pv)
        lleg = lenstr(legend)+1
        if (lleg+lpv-1+2 .gt. maxleg) go to 100
        legend(lleg:lleg+lpv-1+2) = ', '//pv(:lpv)
      end do
 100  continue
      if (dolegn) then
        write (ncsv, '(A)') legend(:lenstr(legend))
      end if

      IF (TIMPLT) THEN
        IF (NPTS(1) .EQ. MAXPTS) THEN
          NX = NTPVAR+1
        ELSE
          NX = NTPVAR+2
        END IF
      ELSE
        NX = 1
      END IF

      do itim = 1, npts(1)
        IF (TIMPLT) THEN
          NY = 1
        ELSE
          NY = 2
        END IF
        write (ncsv, 10100) pts(itim, nx),
     *    (pts(itim,nn),nn=ny,ny+ncrv-1)
10100   FORMAT (1pe15.7E3,100(', ',1pe15.7E3))
      end do
      if (.not. dolegn) then
        write (*,*) 'The data file contains these data:'
        write (*, '(A)') legend(:lenstr(legend))
      else
        write (*,*) NCRV, ' Curves written to CSV file.'
      end if
      write (ncsv,*) ' '
 170  continue
      RETURN
      END

