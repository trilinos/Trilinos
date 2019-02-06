C Copyright (c) 2007-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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
C
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
C

C=======================================================================
      SUBROUTINE DBIVTT (NDB, ISEVOK, ITMP, NELBLK, NVAREL)
C=======================================================================
C
C   --*** DBIVTT *** Read element variable truth table
C   --   Modified for ExodusII format 8/26/95
C   --
C   --Parameters:
C   --   NDB    - IN  - the database number
C   --   NELBLK - IN  - the number of element blocks
C   --   NVAREL - IN  - the number of element variables; <0 if end-of-file
C   --   ISEVOK - OUT - the dynamic memory array of the element block variable
C   --                  truth table;variable i,block j exists iff ISEVOK(j,i)

      INTEGER NDB
      INTEGER NELBLK, NVAREL
      LOGICAL ISEVOK(NELBLK,*)
      INTEGER ITMP(NVAREL, NELBLK)

C     Read the element block variable truth table
C       isevok - num_elem_var cycles faster
      if (nvarel .gt. 0) then
        CALL EXGVTT(NDB, NELBLK, NVAREL, ITMP, IERR)
        DO 110 I = 1, NVAREL
          DO 100 IELB = 1, NELBLK
            ISEVOK(IELB,I) = (ITMP(I,IELB) .NE. 0)
 100      CONTINUE
 110    CONTINUE
      end if
      RETURN
      END

