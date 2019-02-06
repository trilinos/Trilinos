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
      SUBROUTINE DBIST2 (NDB, NVAREL,  NELBLK, NELBDM, ISEVOK,
     $     VAREL, NUMELB, IVAR, IELB, *)
C=======================================================================
C   --*** DBIST2 *** (EXOLIB) Internal to DBISTE, Read element variables
C   --
C   --DBIST2 reads the database element variables for one time step.
C   --
C   --Parameters:
C   --   NDB - IN - the database number
C   --   NVAREL - IN - the number of element variables
C   --   NVARDM - IN - the row dimension of VAREL
C   --   NELBLK - IN - the number of element blocks
C   --   NELBDM - IN - the row dimension of ISEVOK
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   VAREL - OUT - the element variables for the time step (if OPTION)
C   --   IVAR  - OUT - the nodal variable number if error on read.
C   --   IELB  - OUT - the element block number if error on read.
C   --   * - return statement if error encountered, including end-of-file;
C   --      message is printed
C   --
      INTEGER NDB
      INTEGER NVAREL, NELBLK, NELBDM
      INTEGER NUMELB(*)
c      LOGICAL ISEVOK(nvarel,*)
      integer ISEVOK(nvarel,*)
      REAL VAREL(*)

      IELo = 0
      DO 130 IELB = 1, NELBLK
         DO 120 IVAR = 1, NVAREL
            IF (ISEVOK(IVAR,IELB) .ne. 0) THEN
               READ (NDB, END=200, ERR=200, IOSTAT=IERR)
     &              (VAREL(IELo+N), N=1,NUMELB(IELB))
               ielo=ielo+numelb(ielb)
            END IF
  120    CONTINUE
  130 CONTINUE
      RETURN
  200 CONTINUE
      RETURN 1
      END
