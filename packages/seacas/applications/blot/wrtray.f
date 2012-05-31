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

C $Log: wrtray.f,v $
C Revision 1.3  2009/03/25 12:36:49  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  1994/07/21 15:28:20  gdsjaar
C Moved more commons into includes.
C
c Revision 1.1  1994/04/07  20:17:45  gdsjaar
c Initial checkin of ACCESS/graphics/blotII2
c
c Revision 1.3  1993/09/20  18:15:40  gdsjaar
c Changed to output block id instead of number, changed output format slightly
c
c Revision 1.2  1993/09/16  21:13:51  gdsjaar
c Redid method of writing rayshade file. Now, hidden 6 only writes file,
c but does not do shaded plot.  Also, new file written each time rather
c than appending on to end of first file.
c
c
C=======================================================================
      SUBROUTINE wrtray (LENF, NLNKF, LINKF, NXFAC, IXFAC,
     &   XN, YN, ZN, IELBST, BLKCOL, IDELB, *)
C=======================================================================

C   --*** WRTRAY *** (DETOUR) Write polygons to rayshade input file
C   --
C   --WRTRAY writes the polygons to a rayshade input file format.
C   --
C   --Parameters:
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   NXFAC - IN - the number of ordered faces (if DOIXF)
C   --   IXFAC - IN - the indices of the ordered faces (if DOIXF)
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   IELBST - IN - the element block status (>0 if selected)
C   --   BLKCOL - IN/OUT - the user selected colors of the element blocks.
C   --                    BLKCOL(0) = 1 if the user defined material
C   --                                colors should be used in mesh plots.
C   --                              = -1 if program selected colors should
C   --                                be used.
C   --                    BLKCOL(i) = the user selected color of element
C   --                               block i:
C   --                                  -2 - no color selected by user.
C   --                                  -1 - black
C   --                                   0 - white
C   --                                   1 - red
C   --                                   2 - green
C   --                                   3 - yellow
C   --                                   4 - blue
C   --                                   5 - cyan
C   --                                   6 - magenta
C   --   * - return statement if the cancel function is active
C   --
C   --Common Variables:
C   --   Uses NUMEL, NELBLK of /DBNUMS/
C   --   Uses IS3DIM of /D3NUMS/

      include 'dbname.blk'
      include 'dbnums.blk'
      include 'd3nums.blk'

      INTEGER LENF(0:NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      INTEGER IXFAC(*)
      REAL XN(*), YN(*), ZN(*)
      INTEGER IELBST(NELBLK)
      INTEGER BLKCOL(0:NELBLK)
      INTEGER IDELB(*)

      REAL XPTS(0:20), YPTS(0:20), ZPTS(0:20)
      CHARACTER*20 STRING
      CHARACTER*2048 FILNAM
      SAVE NFIL
      DATA NFIL   /0/

C ... Open the output file.  
      filnam = basenam(:lenstr(basenam)) // '.ray'
      NFIL = NFIL + 1
      IF (NFIL .GT. 1) then
        call intstr (1, -1, NFIL, STRING, LSTR)
        FILNAM(LNAM+1:LNAM+LSTR) = STRING(:LSTR)
        LNAM = LNAM + LSTR
      END IF
      open (unit=95, file=filnam(:lnam), form='formatted',
     *  status='unknown')
      
      PI = ATAN2(0.0, -1.0)

      xmin = 1.0e38
      ymin = 1.0e38
      zmin = 1.0e38
      xmax =-1.0e38
      ymax =-1.0e38
      zmax =-1.0e38
      
      LSTBLK = -999
      DO 100 IX = 1, NXFAC
         IFAC = IXFAC(IX)
         IELB = 0
         IXL = IDBLNK (IELB, IFAC, LENF, NLNKF)

         if (lstblk .ne. ielb) then
           write (string, 900) idelb(ielb)
           call pckstr (1, string)
           lstr = lenstr(string)
           lstblk = ielb
         end if
         
         NNPF = NLNKF(IELB)
         do 90 ilink = 0, nnpf-1
           XPTS(ILINK) = XN(LINKF(IXL+ILINK))
           YPTS(ILINK) = YN(LINKF(IXL+ILINK))
           ZPTS(ILINK) = ZN(LINKF(IXL+ILINK))

           xmin = min(xmin, xpts(ilink))
           ymin = min(ymin, ypts(ilink))
           zmin = min(zmin, zpts(ilink))

           xmax = max(xmax, xpts(ilink))
           ymax = max(ymax, ypts(ilink))
           zmax = max(zmax, zpts(ilink))
 90      continue

      write (95, 910) string(:lstr),
     *  (xpts(i), ypts(i), zpts(i),i=0, nnpf-1)

  100 CONTINUE


      write (*,*) 'RAYSHADE polygon file written to ',filnam(:lnam)
      write (*,*)
      write (*,920) 'X', xmin, xmax, xmax-xmin, (xmin+xmax)/2
      write (*,920) 'Y', ymin, ymax, ymax-ymin, (ymin+ymax)/2
      write (*,920) 'Z', zmin, zmax, zmax-zmin, (zmin+zmax)/2
      write (*,*)

      dmax  = max(xmax-xmin, ymax-ymin) / 2.0
      do 110 i=1, 9
        fov = i * 5.0
        angle = FOV/2.0 * PI/180.
        zdist = dmax / tan(angle)
        eyep  = zmax + zdist
        write (*,930) eyep, FOV
 110  continue
        write (*,*)
        close (95)

      RETURN

 900  FORMAT ('id',i10)
 910  format ('poly ',A, 20(1pe13.5))
 920  format (A1,': Min =',1pe11.3,' Max =',1pe11.3,' Range =',
     *  1pe11.3,' Center =',1pe11.3)
 930  FORMAT ('Z-coordinate of eye point =',1pe11.3,' for a ',
     *  0pf4.1,' degree field of view')
      END
