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
      SUBROUTINE RDSPLN (STORE, NDB, SLTOP, SLBOT, NSPL,
     *     ZS, XS, YS)
C=======================================================================

      LOGICAL STORE, MATSTR
      REAL    ZS(*), XS(*), YS(*)
      INTEGER NSPL
      REAL    SLTOP(2), SLBOT(2)
      INTEGER NDB

      PARAMETER (BINGO = 1.0E38)

      PARAMETER (MXFLD = 4)
      REAL RVAL(MXFLD)
      INTEGER KVAL(MXFLD), IVAL(MXFLD)
      CHARACTER*8  CVAL(MXFLD)

      REWIND (NDB)
      IPTA = 0
      SLTOP(1) = BINGO
      SLBOT(1) = BINGO
      SLTOP(2) = BINGO
      SLBOT(2) = BINGO

   10 CONTINUE
      CALL FREFLD ( NDB, 0, 'AUTO', MXFLD, IERR,
     *     NFLD, KVAL, CVAL, IVAL, RVAL)
      IF (IERR .EQ. 0) THEN
         IF (KVAL(1) .EQ. 0) THEN
            IF (MATSTR(CVAL(1), 'TOP', 1) .OR.
     &           MATSTR(CVAL(1), 'FRONT', 1)) THEN
            ELSE IF (MATSTR(CVAL(1), 'BOTTOM', 1) .OR.
     &              MATSTR(CVAL(1), 'BACK', 1)) THEN
            ELSE IF (MATSTR(CVAL(1), 'SLOPE', 1)) THEN
               IF (MATSTR(CVAL(2), 'TOP', 1) .OR.
     &              MATSTR(CVAL(2), 'FRONT', 1)) THEN
                  IF (STORE) THEN
                     SLTOP(1) = RVAL(3)
                     SLTOP(2) = RVAL(4)
                  END IF
               ELSE IF (MATSTR(CVAL(2), 'BACK', 1) .OR.
     &                 MATSTR(CVAL(2), 'BOTTOM', 1)) THEN
                  IF (STORE) THEN
                     SLBOT(1) = RVAL(3)
                     SLBOT(2) = RVAL(4)
                  END IF
               END IF
            END IF
         ELSE IF (NFLD .GE. 2) THEN
            IPTA = IPTA + 1
            IF (STORE) THEN
               ZS(IPTA)   = RVAL(1)
               XS(IPTA)   = RVAL(2)
               YS(IPTA)   = RVAL(3)
            END IF
         END IF
         GO TO 10
      END IF

      IF (.NOT. STORE) THEN
         NSPL = IPTA
      END IF

      RETURN
      END
