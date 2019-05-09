C Copyright(C) 2009-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
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
C     * Neither the name of NTESS nor the names of its
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

C-----------------------------------------------------------------------
      subroutine setcol(nelblk, shdcol, ishdcl, ielbst, blkcol, idelb)
C-----------------------------------------------------------------------
      external blkdat

C ... SHDCOL(1, *) = Red   Component
C     SHDCOL(2, *) = Green Component
C     SHDCOL(3, *) = Blue  Component
C     SHDCOL(4-7,*) - Future Use
      REAL SHDCOL(7,NELBLK)
C ... ISHDCL(1, *) = -1 if color not set, >0 if color set
C     ISHDCL(2, *) = Number of colors to use for this block (SET if 0)
C     ISHDCL(3, *) = Starting location in color map (SET)
      INTEGER ISHDCL(3,NELBLK)
      INTEGER IELBST(NELBLK)
      INTEGER BLKCOL(0:NELBLK)
      INTEGER IDELB(*)
      LOGICAL BLKCLR
      integer IPARMS(3)
      character*80 string

      include 'icrnbw.blk'
      include 'light.blk'

      logical first
      save    first
      data    first /.TRUE./


      if (FIRST) THEN
        FIRST = .FALSE.
        call grcolu('ALTERNATE')
        call grspar('SPECTRUM', -1, 256, STRING)
      end if

      blkclr = .FALSE.
      do 10 i = 1, NELBLK
        if (ielbst(i) .gt. 0 .and. ishdcl(1,i) .gt. 0) then
          blkclr =.TRUE.
          ishdcl(1,i) = i
        end if
 10   continue

      if (blkclr) then

C ... Determine whether any of the block colors match.
        do 18 iblk = 1, nelblk
          if (ishdcl(1,iblk) .eq. iblk .and. ielbst(iblk) .gt. 0) then
            ired = int(10000 * shdcol(1,iblk))
            igrn = int(10000 * shdcol(2,iblk))
            iblu = int(10000 * shdcol(3,iblk))
            do 16 i = iblk+1, nelblk
              if (ielbst(i) .gt. 0) then
                if (ired .eq. int(10000 * shdcol(1,i)) .and.
     *              igrn .eq. int(10000 * shdcol(2,i)) .and.
     *              iblu .eq. int(10000 * shdcol(3,i))) then
                  ishdcl(1,i) = iblk
                  ishdcl(2,iblk) = max(ishdcl(2,iblk), ishdcl(2,i))
                  ishdcl(2,i)    = ishdcl(2,iblk)
                end if
              end if
 16         continue
          end if
 18     continue

        ncol = 0
        icol = 0
        do 20 iblk = 1, NELBLK
          if (ielbst(iblk) .gt. 0 .and. ishdcl(1,iblk) .eq. iblk) then
            if (ishdcl(2, iblk) .eq. 0) then
              icol = icol + 1
            else
              ncol = ncol + ishdcl(2, iblk)
            end if
          end if
 20     continue
C ... If number of colors not set for at least one block,
C       determine number of colors specified vs. spectrum colors,
C       divide remainder of colors to unspecified blocks
        call grgpar('SPECTRUM', 0, IPARMS, STRING)
        NUMCOL = IPARMS(2)
        IF (NCOL .GT. NUMCOL) then
          call prterr ('FATAL','Maximum number of colors exceeded')
        end if
        if (icol .gt. 0) then
          icol = (numcol - ncol) / icol
          ncol = numcol
        else
C ... If number of colors set for all blocks, allocate that many colors.
          call grspar('SPECTRUM', -1, min(NCOL,NUMCOL), STRING)
        end if
        mincol = 0
        do 30 iblk = 1, NELBLK
          if (ielbst(iblk) .gt. 0) then
            if (ishdcl(2, iblk) .eq. 0) ishdcl(2,iblk) = icol
            if (ishdcl(1, iblk) .eq. iblk) then
              ishdcl(3,iblk) = mincol
              mincol = mincol + ishdcl(2, iblk)
            else
              ishdcl(3,iblk) = ishdcl(3,ishdcl(1,iblk))
            end if
          end if
 30     continue

C ... Set the spectrum
        do 50 iblk = 1, nelblk
          if (ielbst(iblk) .gt. 0 .and. ishdcl(1,iblk) .eq. iblk) then
            ncol = ishdcl(2, iblk)
            mincol = ishdcl(3, iblk)
            if (ncol .gt. 1) then
              R = (1.0 / ambient) ** (1.0 / (ncol-1))
            else
              R = 1.0
            end if
            do 40 i = 0, ncol-1
              if (ncol .eq. 1) then
                hue = 1.0
              else
                hue = ambient * R**i
              end if
              red   = min(1.0, hue * shdcol(1, iblk) )
              green = min(1.0, hue * shdcol(2, iblk) )
              blue  = min(1.0, hue * shdcol(3, iblk) )
              CALL PLTCOL (i+mincol+8, red, green, blue)
 40         continue
          end if
 50     continue
        call prnshd (nelblk, idelb, ishdcl, shdcol, ielbst)
      else
C ... Single color for entire model
        mincol = 0
        call grgpar('SPECTRUM', 0, IPARMS, STRING)
        NCOL = IPARMS(1)
        do 80 iblk = 1, nelblk
          ishdcl(2,iblk) = ncol
          ishdcl(3,iblk) = mincol
 80     continue
        if (ncol .gt. 1) then
          R = (1.0 / ambient) ** (1.0 / (ncol-1))
        else
          R = 1.0
        end if
        do 90 i=0, NCOL-1
          hue = ambient * R**i
          red   = min(1.0, (hue * RMULT))
          green = min(1.0, (hue * GMULT))
          blue  = min(1.0, (hue * BMULT))
          CALL PLTCOL (i+8, red, green, blue)
 90     continue
      end if
      call grcolu('ALTERNATE')

      return
      end
C-----------------------------------------------------------------------
C $Log: setcol.f,v $
C Revision 1.4  2009/03/25 12:36:47  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.3  1998/08/21 14:57:59  gdsjaar
C Fixed array bounds problem
C
C Revision 1.2  1998/06/12 15:53:36  gdsjaar
C 1. Problem with TIMES array. Blot accesses a dummy timestep even if
C there were no timesteps on the database. Array wasn't allocated, so
C writing off into never-never land.
C
C 2. Inconsistency among some common blocks. Some places weren't using
C the include but had the definition hardwired in. Removed those.
C
C 3. Added 'EXTERNAL BLKDAT' to all routines that used data values set
C in BLKDAT
C
C 4. Cleanup of some A vs. IA argument passing.
C
C Revision 1.1  1994/04/07 20:11:29  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.1  1993/10/22  18:56:42  gdsjaar
c Moved block shading routine to separate subroutine.
c

