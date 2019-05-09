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
      SUBROUTINE NEWATR (NELBLK, NUMATR, ATRSC, NUMELB, ATRIB)
C=======================================================================

      INTEGER NUMATR(*)
      INTEGER NUMELB(*)
      REAL ATRSC(2,*)
      REAL ATRIB(*)

      IEATR = 0
      IAT   = 1

      DO 100 IELB = 1, NELBLK
         ISATR = IEATR + 1
         IEATR = IEATR + NUMATR(IELB) * NUMELB(IELB)
         if (numatr(ielb) .gt. 0) then
           CALL NEWAT1 (NUMELB(IELB), NUMATR(IELB),
     *       ATRSC(1,IAT), ATRIB(ISATR))
           IAT = IAT + NUMATR(IELB)
         end if
  100 CONTINUE

      RETURN
      END

      SUBROUTINE NEWAT1(NUMEL, NUMATR, ATRSC, ATRIB)
      REAL ATRSC(2,*)
      REAL ATRIB(*)

      IBEG = 1
      DO 110 IATR = 1, NUMATR
        if (ATRSC(1,IATR) .NE. 0.0) then
          DO 100 IEL = 1, NUMEL
            ATRIB(IBEG+NUMATR*(IEL-1)) = ATRSC(1,IATR)
 100      CONTINUE
        else if (ATRSC(2,IATR) .NE. 1.0) then
          DO 105 IEL = 1, NUMEL
            ATRIB(IBEG+NUMATR*(IEL-1)) = ATRSC(2,IATR) *
     *        ATRIB(IBEG+NUMATR*(IEL-1))
 105      CONTINUE
        end if
        IBEG = IBEG + 1
 110  CONTINUE

      RETURN
      END

