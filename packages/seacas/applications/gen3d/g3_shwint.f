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

C=======================================================================
      SUBROUTINE SHWINT (ITRANT, NEREPL, DIM3, NRTRAN, D3TRAN, ZGRAD)
C=======================================================================

      INTEGER NRTRAN(*)
      REAL    D3TRAN(*), ZGRAD(*)

      CHARACTER*20 RSTR(9)
      CHARACTER*20 STRA, TYPE

      CALL INTSTR (1, 0, NEREPL, STRA, LSTRA)
      CALL NUMSTR (1, 4, DIM3, RSTR(1), LR)

      IF      (ITRANT .EQ.  0) THEN
         TYPE = 'Transform'
      ELSE IF (ITRANT .EQ.  1) THEN
         TYPE = 'Translate'
      ELSE IF (ITRANT .EQ.  2) THEN
         TYPE = 'Rotate'
      ELSE IF (ITRANT .EQ.  4) THEN
         TYPE = 'Warp'
      ELSE IF (ITRANT .EQ.  8) THEN
         TYPE = 'Twist'
      ELSE IF (ITRANT .EQ. 16) THEN
         TYPE = 'Project'
      ELSE IF (ITRANT .EQ. 32) THEN
         TYPE = 'ExpRotate'
      ELSE IF (ITRANT .EQ. 64) THEN
         TYPE = 'Spline'
      ELSE
         CALL PRTERR ('PROGRAM', 'Unknown transformation option')
         RETURN
      END IF

      LT = LENSTR(TYPE)

      IF (NEREPL .EQ. NRTRAN(1)) THEN
         IF (ABS (ZGRAD(1) - 1.0) .LE. 1.0E-6) THEN
            WRITE (*, 20) TYPE(:LT), ' mesh ', STRA(:LSTRA),
     &         ' times for a total of ', RSTR(1)(:LR)
         ELSE
            CALL NUMSTR (1, 3, ZGRAD(1), RSTR(2), LR2)
            WRITE (*, 20) TYPE(:LT), ' mesh ', STRA(:LSTRA),
     &         ' times for a total of ', RSTR(1)(:LR),
     &         ' with a gradient of ', RSTR(2)(:LR2)
         END IF
      ELSE
         IBLK = 0
         NR = 0
   10    CONTINUE
         IF (.TRUE.) THEN
            IBLK = IBLK + 1
            NR = NR + NRTRAN(IBLK)
            CALL INTSTR (1, 0, NRTRAN(IBLK), STRA, LSTRA)
            CALL NUMSTR (1, 4, D3TRAN(IBLK), RSTR(1), LR)
            IF (ABS (ZGRAD(IBLK) - 1.0) .LE. 0.001) THEN
               WRITE (*, 20) TYPE(:LT), ' mesh ',
     &            STRA(:LSTRA), ' times for a subtotal of ',
     &            RSTR(1)(:LR)
            ELSE
               CALL NUMSTR (1, 3, ZGRAD(IBLK), RSTR(2), LR2)
               WRITE (*, 20) TYPE(:LT), ' mesh ',
     &            STRA(:LSTRA), ' times for a subtotal of ',
     &            RSTR(1)(:LR), ' with a gradient of ',
     &            RSTR(2)(:LR2)
            END IF
            IF (NR .LT. NEREPL) GOTO 10
         END IF
      END IF
   20 FORMAT (1X, 10A)
      RETURN
      END
