C    Copyright(C) 2011 Sandia Corporation.  Under the terms of Contract
C    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C    certain rights in this software
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C              
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C                            
C    * Neither the name of Sandia Corporation nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
C                                                    
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C    
C=======================================================================
      SUBROUTINE PRFRM (NOUT)
C=======================================================================

C   --*** PRFRM *** (GROPE) Display Coordinate Frame information
C   --
C   --Parameters:
C   --   NOUT - IN - the output file, <=0 for standard

      include 'dbase.blk'
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
      include 'dbase.blk'
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
      WRITE(FMT, 35) 'I10', FMT1(:LFMT), FMT1(:LFMT), FMT1(:LFMT)

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
