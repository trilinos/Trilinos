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

c
C=======================================================================
      SUBROUTINE PRNAME (NOUT, NAMLEN, 
     *                   NVARGL, NVARNP, NVAREL, NVARNS, NVARSS,
     &                   NAMEGV, NAMENV, NAMEEV, NAMNSV, NAMSSV)
C=======================================================================

C   --*** PRNAME *** (BLOT) Display database variable names
C   --   Written by Amy Gilkey - revised 01/14/88
C   --
C   --PRNAME displays the database variable names.
C   --
C   --Parameters:
C   --   NOUT - IN - the output file, <=0 for standard
C   --   NVARGL - IN - the number of global variables
C   --   NVARNP - IN - the number of nodal variables
C   --   NVAREL - IN - the number of element variables
C   --   NAMEGV - IN - the global variable names
C   --   NAMENV - IN - the nodal variable names
C   --   NAMEEV - IN - the element variable names

      CHARACTER*(NAMLEN) NAMEGV(*)
      CHARACTER*(NAMLEN) NAMENV(*)
      CHARACTER*(NAMLEN) NAMEEV(*)
      CHARACTER*(NAMLEN) NAMNSV(*)
      CHARACTER*(NAMLEN) NAMSSV(*)
      CHARACTER*128 FMT1, FMT

      IF (NOUT .GT. 0) WRITE (NOUT, 10000)

      IF (NOUT .GT. 0) THEN
         WRITE (NOUT, 10010)
      ELSE
         WRITE (*, 10010)
      END IF

      WRITE(FMT1,20) NAMLEN
      CALL SQZSTR(FMT1, LFMT)
      if (namlen .le. 20) then
         WRITE(FMT, 30) FMT1(:LFMT), FMT1(:LFMT)
      else
         WRITE(FMT, 40) FMT1(:LFMT), FMT1(:LFMT)
      endif

C ... Print them out.
      if (nout .le. 0) then
            WRITE (*, FMT) 'Global: ', (NAMEGV(I), I=1,NVARGL)
            WRITE (*, FMT) 'Nodal:  ', (NAMENV(I), I=1,NVARNP)
            WRITE (*, FMT) 'Element:', (NAMEEV(I), I=1,NVAREL)
            WRITE (*, FMT) 'Nodeset:', (NAMNSV(I), I=1,NVARNS)
            WRITE (*, FMT) 'Sideset:', (NAMSSV(I), I=1,NVARSS)
      else
            WRITE (NOUT, FMT) 'Global: ', (NAMEGV(I), I=1,NVARGL)
            WRITE (NOUT, FMT) 'Nodal:  ', (NAMENV(I), I=1,NVARNP)
            WRITE (NOUT, FMT) 'Element:', (NAMEEV(I), I=1,NVAREL)
            WRITE (NOUT, FMT) 'Nodeset:', (NAMNSV(I), I=1,NVARNS)
            WRITE (NOUT, FMT) 'Sideset:', (NAMSSV(I), I=1,NVARSS)
      end if

      RETURN

 20   FORMAT('A',I4)
 30   FORMAT ('(4X, A, :, 3 (2X, ',A,'), :, /,(12X, 3 (2X, ',A,')))')
 40   FORMAT ('(4X, A, :, 2 (2X, ',A,'), :, /,(12X, 2 (2X, ',A,')))')

10000  FORMAT (/, 1X, 'VARIABLES NAMES')
10010  FORMAT (/, 1X, 'Variables Names:')

      END
