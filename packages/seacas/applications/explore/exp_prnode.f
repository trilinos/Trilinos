C    Copyright(C) 2008-2017 National Technology & Engineering Solutions of
C    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
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
C    * Neither the name of NTESS nor the names of its
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
      SUBROUTINE PRNODE (OPTION, NOUT, NUMNP, LISNP,
     &   NVARNP, LISNV, NAMENV, VARNP, MAP, DOMAP)
C=======================================================================

C   --*** PRNODE *** (EXPLORE) Display current database nodal variables
C   --
C   --PRNODE displays the nodal data for a time step.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options
C   --   NOUT - IN - the output file, <=0 for standard
C   --   NUMNP - IN - the number of nodes
C   --   LISNP - IN - the indices of the selected nodes
C   --   NVARNP - IN - the number of nodal variables
C   --   LISNV - IN - the indices of the selected nodal variables
C   --   NAMENV - IN - the names of the nodal variables
C   --   VARNP - IN - the nodal variables for the time step

      include 'exodusII.inc'
      CHARACTER*(*) OPTION
      INTEGER LISNP(0:*)
      INTEGER LISNV(0:*)
      CHARACTER*(*) NAMENV(*)
      REAL VARNP(numnp,*)
      INTEGER MAP(*)
      LOGICAL DOMAP
      INTEGER GETPRC, PRTLEN
      CHARACTER*128 FMT1, FMT

      PRTLEN = GETPRC() + 7
      WRITE(FMT1,20) PRTLEN, PRTLEN-7
      CALL SQZSTR(FMT1, LFMT)
      WRITE(FMT, 30) FMT1(:LFMT), FMT1(:LFMT)

      IF (NOUT .GT. 0) WRITE (NOUT, 10000)

      if (domap) then
        if (nout .gt. 0) then
          write (nout, 10005)
        else
          write (*, 10005)
        end if
      end if

      do 90 i=1, lisnv(0)
        irow = ((i-1)/5)+1
        icol = i - (irow-1)*5
        IF (NOUT .GT. 0) THEN
           WRITE (NOUT, 10010) irow, icol, NAMENV(LISNV(I))
        ELSE
           WRITE (*, 10010) irow, icol, NAMENV(LISNV(I))
        END IF
 90   continue


      DO 100 IX = 1, LISNP(0)
         INP = LISNP(IX)
         if (domap) then
            id = map(inp)
         else
            id = inp
         end if

         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, FMT, IOSTAT=IDUM)
     &           ID, (VARNP(INP, LISNV(I)), I=1,LISNV(0))
         ELSE
            WRITE (*, FMT, IOSTAT=IDUM)
     &           ID, (VARNP(INP, LISNV(I)), I=1,LISNV(0))
         END IF
 100  CONTINUE

      RETURN

 20   FORMAT('1PE',I2.2,'.',I2.2)
 30   FORMAT('(1X, ''Node'', I12, 5 (2X, ',A,'), :, /,',
     $     '(15X, 5 (2X, ',A,')))')
10000 FORMAT (/, 1X, 'NODAL TIME STEP VARIABLES')
10005 FORMAT (1X, 'Nodal ids are Global')
10010 FORMAT (1X, 'Row ',I4,', Column ',I1,' is variable ',A)
      END
