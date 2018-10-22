C Copyright(C) 2011-2017 National Technology & Engineering Solutions of
C Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C           
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C                         
C * Neither the name of NTESS nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
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
      SUBROUTINE WRPXYZ (X, Y, Z, NUMNP, IWARP, NORMAL, REFDIS)
C=======================================================================
      REAL X(NUMNP), Y(NUMNP), Z(NUMNP)

C ... Origin
      if (iwarp .eq. 1) then
        call PRTERR('PROGRAM', 'Origin warping not implemented yet.')
C ... Xaxis
      else if (iwarp .eq. -1) then
        if (normal .eq. 2) then
          call warpit(x, y, z, numnp, refdis)
        else if (normal .eq. 3) then
          call warpit(x, z, y, numnp, refdis)
        end if
        
C ... Yaxis
      else if (iwarp .eq. -2) then
        if (normal .eq. 3) then
          call warpit(y, z, x, numnp, refdis)
        else if (normal .eq. 1) then
          call warpit(y, x, z, numnp, refdis)
        end if
        
C ... Zaxis
      else if (iwarp .eq. -3) then
        if (normal .eq. 1) then
          call warpit(z, x, y, numnp, refdis)
        else if (normal .eq. 2) then
          call warpit(z, y, x, numnp, refdis)
        end if
      end if
      
      RETURN
      END
      
      SUBROUTINE WARPIT(C1, C2, C3, NUMNP, REFDIS)
      REAL C1(NUMNP), C2(NUMNP), C3(NUMNP)
      
      do 10 i=1, numnp
        c1(i) = c1(i)
        
        radius = c2(i)
        theta  = c3(i) / refdis
        
        c3(i) = radius * sin(theta)
        c2(i) = radius * cos(theta)
 10   continue
      
      return
      end
