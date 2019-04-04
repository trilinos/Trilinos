C Copyright(C) 2011-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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

C -*- Mode: fortran -*-
      SUBROUTINE RDSPLN (DOUBLE, STORE, NDB, RDTHET, SLLFT, SLRGT, NSPL,
     *   RSA, ZSA, RSB, ZSB)

C ... If only single spline, then RSB is Thickness

      LOGICAL DOUBLE, STORE
      REAL    RSA(*), ZSA(*), RSB(*), ZSB(*)
      INTEGER NSPL(2)
      REAL    SLLFT(2), SLRGT(2)
      INTEGER NDB
      LOGICAL RDTHET, MATSTR

      PARAMETER (BINGO = 1.0E38)

      PARAMETER (MXFLD = 4)
      REAL RVAL(MXFLD)
      INTEGER KVAL(MXFLD), IVAL(MXFLD)
      CHARACTER*8  CVAL(MXFLD)

      REWIND (NDB)
      RDTHET = .FALSE.
      IPTA = 0
      IPTB = 0
      SLLFT(1) = BINGO
      SLRGT(1) = BINGO
      SLLFT(2) = BINGO
      SLRGT(2) = BINGO
      ISPL = 1

   10 CONTINUE
      CALL FREFLD ( NDB, 0, 'AUTO', MXFLD, IERR,
     *   NFLD, KVAL, CVAL, IVAL, RVAL)
      IF (IERR .EQ. 0) THEN
         IF (KVAL(1) .EQ. 0) THEN
            IF (MATSTR(CVAL(1), 'LEFT', 1)) THEN
               if (store) SLLFT(ISPL) = RVAL(2)
            ELSE IF (MATSTR(CVAL(1), 'RIGHT', 1)) THEN
               if (store) SLRGT(ISPL) = RVAL(2)
            ELSE IF (MATSTR(CVAL(1), 'ANGULAR', 1)) THEN
               RDTHET = .TRUE.
            ELSE IF (MATSTR(CVAL(1), 'TOP', 1) .OR.
     &         MATSTR(CVAL(1), 'FRONT', 1)) THEN
               IF (.NOT. STORE) DOUBLE = .TRUE.
               ISPL = 1
            ELSE IF (MATSTR(CVAL(1), 'BOTTOM', 1) .OR.
     &         MATSTR(CVAL(1), 'BACK', 1)) THEN
               IF (.NOT. STORE) DOUBLE = .TRUE.
               ISPL = 2
            else if (matstr(cval(1), 'slope', 1)) then
               if (matstr(cval(2), 'top', 1) .or.
     &             matstr(cval(2), 'front', 1)) then
                   if (store) then
                      sllft(1) = rval(3)
                      sllft(2) = rval(4)
                   end if
               else if (matstr(cval(2), 'back', 1) .or.
     &             matstr(cval(2), 'bottom', 1)) then
                   if (store) then
                      slrgt(1) = rval(3)
                      slrgt(2) = rval(4)
                   end if
               end if
            END IF
         ELSE IF (NFLD .GE. 2) THEN
            IF (ISPL .EQ. 1) THEN
               IPTA = IPTA + 1
               IF (STORE) THEN
                  RSA(IPTA)   = RVAL(1)
                  ZSA(IPTA)   = RVAL(2)
                  IF (.NOT. DOUBLE) RSB(IPTA) = RVAL(3)
               END IF
            ELSE
               IPTB = IPTB + 1
               IF (STORE) THEN
                  RSB(IPTB)   = RVAL(1)
                  ZSB(IPTB)   = RVAL(2)
               END IF
            END IF
         END IF
         GO TO 10
      END IF

      IF (.NOT. STORE) THEN
         NSPL(1) = IPTA
         NSPL(2) = IPTB
      END IF

      RETURN
      END
