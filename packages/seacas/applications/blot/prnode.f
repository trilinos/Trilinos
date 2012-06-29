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

C=======================================================================
      SUBROUTINE PRNODE (OPTION, NOUT, NUMNP, LISNP,
     &   NVARNP, LISNV, NAMENV, VARNP, MAPND)
C=======================================================================

C   --*** PRNODE *** (BLOT) Display current database nodal variables
C   --   Written by Amy Gilkey - revised 11/05/87
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
C   --   VARNP - IN - the selected nodal variables for the time step

      CHARACTER*(*) OPTION
      INTEGER LISNP(0:*)
      INTEGER LISNV(0:*)
      CHARACTER*(*) NAMENV(*)
      REAL VARNP(NUMNP,*)
      INTEGER MAPND(*)
      
      LOGICAL ISABRT

      IF (NOUT .GT. 0) WRITE (NOUT, 10000)

      do i=1, lisnv(0)
        irow = ((i-1)/5)+1
        icol = i - (irow-1)*5
        IF (NOUT .GT. 0) THEN
          WRITE (NOUT, 10010) irow, icol, NAMENV(LISNV(I))
        ELSE
           WRITE (*, 10010) irow, icol, NAMENV(LISNV(I))
        END IF
      end do

      DO 100 IX = 1, LISNP(0)
         IF (ISABRT ()) RETURN
         INP = LISNP(IX)
         IF (NOUT .GT. 0) THEN
            WRITE (NOUT, 10020, IOSTAT=IDUM)
     &       MAPND(INP), (VARNP(INP,I), I=1,LISNV(0))
         ELSE
            WRITE (*, 10020, IOSTAT=IDUM)
     &       MAPND(INP), (VARNP(INP,I), I=1,LISNV(0))
         END IF
  100 CONTINUE

      RETURN

10000 FORMAT (/, 1X, 'NODAL TIME STEP VARIABLES')
10010 FORMAT (1X, 'Row ',I4,', Column ',I1,' is variable ',A)
10020 FORMAT (1X, 'Node', I6, 5 (1X, 1PE13.6), :, /,
     &   (1X, 10X, 5 (1X, 1PE13.6)))
      END
