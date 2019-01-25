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
      subroutine prnshd (nelblk, idelb, ishdcl, shdcol, ielbst)
C-----------------------------------------------------------------------
      INTEGER IDELB(*)
      INTEGER ISHDCL(3,NELBLK)
      REAL SHDCOL(7,NELBLK)
      INTEGER IELBST(NELBLK)

      CHARACTER*132 ELBLIN

      do 110 iblk = 1, nelblk
        if (ishdcl(1, iblk) .eq. iblk) then
          if (ishdcl(2, iblk) .gt. 1) then
            elblin = ' Shades'
          else
            elblin = ' Shade'
          end if
          write (*, 920) (shdcol(i,iblk),i=1,6), ishdcl(2, iblk),
     *      elblin(:lenstr(elblin))
          ELBLIN = ' Block IDs: '
          N = 2
          do 100 i = iblk, nelblk
            if (ishdcl(1, i) .eq. ishdcl(1, iblk)) then
              n = n + 1
              WRITE (ELBLIN((N-1)*6+1:N*6), '(I6)', IOSTAT=IDUM)
     *          idelb(i)
              IF (N .GE. 12) THEN
                WRITE (*, 940) ELBLIN(1:LENSTR(ELBLIN))
                N = 0
                ELBLIN = ' '
              END IF
            end if
 100      continue
          IF (N .GT. 0) THEN
            WRITE (*, 940) ELBLIN(1:LENSTR(ELBLIN))
            write (*,*) ' '
          END IF
        end if
 110  continue

 920  FORMAT (' Block Shading: RGB =', 6(f6.3),', ',i3,A)
 940  FORMAT (A)

      RETURN
      end


