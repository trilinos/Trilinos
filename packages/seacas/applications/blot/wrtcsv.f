C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
C         
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
      call flush(ncsv)
 170  continue
      RETURN
      END
      
      
