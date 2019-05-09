C Copyright (C) 2009-2017 National Technology & Engineering Solutions
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
      SUBROUTINE DBONM1 (NDB, NELBLK, NVAREL, ISEVOK, IEVOK, NBLKDM)
C=======================================================================
C   --*** DBONM1 *** (EXOLIB) Internal to DBONAM
C   --
C   --DBONM1 writes the element block variable truth table.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   NELBLK - IN - the number of element blocks
C   --   NVAREL - IN - the number of element variables
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(i,j)
C   --   IEVOK - SCRATCH - size = NELBLK * NVAREL (may overlap ISEVOK)
C   --
C   --Database must be positioned in front of truth table upon entry;
C   --upon exit positioned after table.

      INTEGER NDB
      INTEGER NELBLK, NVAREL
c      LOGICAL ISEVOK(nvarel,*)
      integer ISEVOK(nvarel,*)
      INTEGER IEVOK(nvarel,*)

C   --Write the element block variable truth table

      IF ((NVAREL .GT. 0) .AND. (NELBLK .GT. 0)) THEN
         DO 110 IELB = 1, NELBLK
            DO 100 I = 1, NVAREL
               IF (ISEVOK(I,IELB) .ne. 0) THEN
                  IEVOK(I,IELB) = 1
               ELSE
                  IEVOK(I,IELB) = 0
               END IF
  100       CONTINUE
  110    CONTINUE

         WRITE (NDB) ((IEVOK(I,IELB), I=1,NVAREL), IELB=1,NELBLK)

c         DO 130 IELB = 1, NELBLK
c            DO 120 I = 1, NVAREL
c               ISEVOK(IELB,I) = (IEVOK(IELB,I) .EQ. 1)
c  120       CONTINUE
c  130    CONTINUE
      ELSE
         WRITE (NDB) 0
      END IF

      RETURN
      END
