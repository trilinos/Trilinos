C    Copyright(C) 2009-2017 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C
C        * Redistributions of source code must retain the above copyright
C          notice, this list of conditions and the following disclaimer.
C
C        * Redistributions in binary form must reproduce the above
C          copyright notice, this list of conditions and the following
C          disclaimer in the documentation and/or other materials provided
C          with the distribution.
C
C        * Neither the name of NTESS nor the names of its
C          contributors may be used to endorse or promote products derived
C          from this software without specific prior written permission.
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

      subroutine hunt(a,n,x,jlo)

C ... find jlo such that a(jlo) .le. x .and. a(jlo+1) .gt. x
C     (if jlo .ne. n)
C
C     Start search at passed in 'jlo' position
C
      DIMENSION a(n)

      integer jlo, low, high

      if (jlo .lt. 1 .or. jlo .gt. n) jlo = 1

      if (a(jlo) .eq. x) then
        return
      else if (a(1) .ge. x) then
        jlo = 1
        return
      else if (a(n) .lt. x) then
        jlo = n
        return
      end if

      if (a(jlo) .le. x) then
        low  = jlo
        high = n-1
      else
        low  = 1
        high = jlo
      end if

      do jlo = low, high
        if (a(jlo) .le. x .and. a(jlo+1) .ge. x) then
          return
        end if
      end do
C ... should not get here since extremes checked above...

      end
