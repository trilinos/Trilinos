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

      SUBROUTINE textur(sat, ncolors, map, invert, rmult, gmult, bmult)

      include 'colormap.blk'
      INTEGER invert, map, ncolors, x, ihue
      REAL hue, sat, red, green, blue, h

      if (ncolors .eq. 0) ncolors = 7

      sat = min(sat, 0.99)
C ... Saturation is used to pass in the minimum hue for some maps.
      if (sat .ge. 0.90 .or. sat .le. 0.0) then
        huemin = 0.3
      else
        huemin = sat
      end if

      ratio = (1.0 / huemin) ** (1.0 / (ncolors-1))
      do x = 0, ncolors-1
        hue  = float(x) / float(ncolors)
        ihue = int(ncolors*hue)
C ... Linear hue map
        hue  = float(ihue)/float(ncolors-1)
        hue  = min(hue, 1.0)
C ... Logarithmic hue map
        huel = huemin * ratio**x
        huel = min(huel, 1.0)

        if(invert .eq. 1) then
          hue = 1.-hue
          huel= 1.-huel
        end if

        if (map .eq. RAINBW) THEN
          huetmp = 1.0 - hue
          call rainbow(huetmp, sat, 1., red, green, blue)
        else if (map .eq. VIRDIS) THEN
          call viridis(x, ncolors, red, green, blue)
        else if (map .eq. TERRAIN) THEN
          h = 3*hue
          if(h .LT. 0.25) THEN
            red   = 0.0
            green = 0.0
            blue  = 0.25+2*h
          else if(h .LT. 2) THEN
            red   = 0.0
            green = 0.25+(2-h)
            blue  = 0.0
          else if(h .LT. 2.7) THEN
            red   = 0.75
            green = 0.15
            blue  = 0.0
          else
            red   = 0.9
            green = 0.9
            blue  = 0.9
          end if
        else if (map .eq. IRON) then
          red   = 3*(hue + 0.03)
          green = 3*(hue - 1.0/3.0)
          blue  = 3*(hue - 2.0/3.0)
        else if (map .eq. ASTRO) then
          red   = hue
          green = hue
          blue  = (hue+.2)/1.2
        else if (map .eq. ZEBRA) then
          red   = MOD(ihue+invert, 2)
          green = red
          blue  = red
        else if (map .eq. GRAY) THEN
          red   = huel
          green = huel
          blue  = huel
        else if (map .eq. METAL) then
          red   = huel * RMULT
          green = huel * GMULT
          blue  = huel * BMULT
        else if (map .eq. COOL) then
          red   = hue
          green = 1.0 - hue
          blue  = 1.0
        end if
        red   = max(0.0, min(red,   1.0))
        green = max(0.0, min(green, 1.0))
        blue  = max(0.0, min(blue,  1.0))
        CALL PLTCOL (8+x, red, green, blue)
      end do
      return
      end

