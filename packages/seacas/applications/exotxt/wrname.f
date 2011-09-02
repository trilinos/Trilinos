C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
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
C 
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
C 


C=======================================================================
      SUBROUTINE WRNAME (NTXT, NDIM, NELBLK, NVARGL, NVARNP, 
     &                   NVAREL, NAMEGV, NAMENV, NAMEEV, ISEVOK, NAMLEN)
C=======================================================================
C   --
C   --*** WRNAME *** (TXTEXO) Write database names
C   --   Written by Amy Gilkey - revised 02/08/88
C   --   Modified for ExodusIIv2 format 10/16/95
C   --
C   --WRNAME writes the names of the coordinates, the element block types,
C   --and the database variables.  The element block variable truth table
C   --is also written.
C   --
C   --Parameters:
C   --   NTXT   - IN - the text file
C   --   NDIM   - IN - the number of coordinates per node; written only if >= 0
C   --   NELBLK - IN - the number of element blocks; written only if >= 0
C   --   NVARGL - IN - the number of global variables
C   --   NVARNP - IN - the number of nodal variables
C   --   NVAREL - IN - the number of element variables
C   --   NAMEGV - IN - the global variable names
C   --   NAMENV - IN - the nodal variable names
C   --   NAMEEV - IN - the element variable names
C   --   ISEVOK - IN - the element block variable truth table;
C   --                 variable i of block j exists iff ISEVOK(j,i)

      INTEGER NTXT, NDIM, NELBLK, NVARGL, NVARNP, NVAREL
      CHARACTER*(namlen) NAMEGV(*)
      CHARACTER*(namlen) NAMENV(*)
      CHARACTER*(namlen) NAMEEV(*)
      LOGICAL ISEVOK(NELBLK,*)

      IF (NDIM .LT. 0) GOTO 110

C   --Write variable names

      WRITE (NTXT, '(A)') '! Variable names'
      WRITE (NTXT, 10010) NVARGL, NVARNP, NVAREL,
     &   '! global, nodal, element variables'

      WRITE (NTXT, 10000) (NAMEGV(I), I=1,NVARGL)

      WRITE (NTXT, 10000) (NAMENV(I), I=1,NVARNP)

      WRITE (NTXT, 10000) (NAMEEV(I), I=1,NVAREL)

C   --Write the element block variable truth table

      if (nvarel .gt. 0) then
        WRITE (NTXT, '(A)') '! Element block variable truth table'
        DO 100 IELB = 1, NELBLK
          WRITE (NTXT, 10020) (ISEVOK(IELB,I), I=1,NVAREL)
 100    CONTINUE
      end if
  110 CONTINUE
      RETURN
10000  FORMAT (2 (A, 1X))
10010  FORMAT (3I10, 6X, A)
10020  FORMAT (40L2)
      END
