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
      SUBROUTINE PRNAME (OPTION, NOUT, NAMLEN, 
     *                   NVARHI, NVARGL, NVARNP, NVAREL, NVARNS, NVARSS,
     &                   NAMEHV, NAMEGV, NAMENV, NAMEEV, NAMNSV, NAMSSV)
C=======================================================================

C   --*** PRNAME *** (BLOT) Display database variable names
C   --   Written by Amy Gilkey - revised 01/14/88
C   --
C   --PRNAME displays the database variable names.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options:
C   --                 'H' to print history variable names
C   --                 'G' to print global variable names
C   --                 'N' to print nodal variable names
C   --                 'E' to print element variable names
C   --                 'M' to print nodeset variable names
C   --                 'S' to print sideset variable names
C   --   NOUT - IN - the output file, <=0 for standard
C   --   NVARHI - IN - the number of history variables
C   --   NVARGL - IN - the number of global variables
C   --   NVARNP - IN - the number of nodal variables
C   --   NVAREL - IN - the number of element variables
C   --   NAMEHV - IN - the history variable names
C   --   NAMEGV - IN - the global variable names
C   --   NAMENV - IN - the nodal variable names
C   --   NAMEEV - IN - the element variable names

      CHARACTER*(*) OPTION
      CHARACTER*(NAMLEN) NAMEHV(*)
      CHARACTER*(NAMLEN) NAMEGV(*)
      CHARACTER*(NAMLEN) NAMENV(*)
      CHARACTER*(NAMLEN) NAMEEV(*)
      CHARACTER*(NAMLEN) NAMNSV(*)
      CHARACTER*(NAMLEN) NAMSSV(*)

      IF (NOUT .GT. 0) WRITE (NOUT, 10000)

      IF (NOUT .GT. 0) THEN
         WRITE (NOUT, 10010)
      ELSE
         WRITE (*, 10010)
      END IF

C ... Print them out.
      IF (INDEX (OPTION, 'H') .GT. 0) THEN
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10020) 'History:', (NAMEHV(I), I=1,NVARHI)
         ELSE
            WRITE (*, 10020) 'History:', (NAMEHV(I), I=1,NVARHI)
         END IF
      END IF
      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'G') .GT. 0)) THEN
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10020) 'Global: ', (NAMEGV(I), I=1,NVARGL)
         ELSE
            WRITE (*, 10020) 'Global: ', (NAMEGV(I), I=1,NVARGL)
         END IF
      END IF
      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'N') .GT. 0)) THEN
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10020) 'Nodal:  ', (NAMENV(I), I=1,NVARNP)
         ELSE
            WRITE (*, 10020) 'Nodal:  ', (NAMENV(I), I=1,NVARNP)
         END IF
      END IF
      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'E') .GT. 0)) THEN
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10020) 'Element:', (NAMEEV(I), I=1,NVAREL)
         ELSE
            WRITE (*, 10020) 'Element:', (NAMEEV(I), I=1,NVAREL)
         END IF
      END IF
      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'M') .GT. 0)) THEN
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10020) 'Nodeset:', (NAMNSV(I), I=1,NVARNS)
         ELSE
            WRITE (*, 10020) 'Nodeset:', (NAMNSV(I), I=1,NVARNS)
         END IF
      END IF
      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'S') .GT. 0)) THEN
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10020) 'Sideset:', (NAMSSV(I), I=1,NVARSS)
         ELSE
            WRITE (*, 10020) 'Sideset:', (NAMSSV(I), I=1,NVARSS)
         END IF
      END IF

      RETURN

10000  FORMAT (/, 1X, 'VARIABLES NAMES')
10010  FORMAT (/, 1X, 'Variables Names:')
10020  FORMAT (4X, A, :, 2 (2X, A), :, /,
     &   (12X, 2 (2X, A)))

      END
