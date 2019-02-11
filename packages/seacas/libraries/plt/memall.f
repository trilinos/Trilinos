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

C $Id: memall.f,v 1.1 1993/07/16 16:46:56 gdsjaar Exp $
C $Log: memall.f,v $
C Revision 1.1  1993/07/16 16:46:56  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      INTEGER FUNCTION MEMALL(LENGTH,MEMRY)
      INTEGER LENGTH
      INTEGER MEMRY(*)
      INTEGER LR
      INTEGER LMX
      INTEGER IP
      INTEGER IPN
      INTEGER LAV
      INTEGER LT
      INTEGER IB

      IF (LENGTH.LE.0) THEN
         WRITE (6,*) ' Cannot allocate a segment of length zero.'
         CALL CPUTBK(.FALSE.)
      END IF

      LR = LENGTH + MOD(LENGTH,2)
      LMX = MEMRY(1)
      IB = 0
      IP = 3
 2650 IF (.NOT. (IB.EQ.0)) GO TO 2660
      IPN = MEMRY(IP)
      IF (IPN.LE.0) THEN
         LT = IP + LR + 1
         IF (LT.GE.LMX) THEN
            WRITE (6,*) ' Cannot allocate space in memall.'
            WRITE (6,*) ' lt >= lmx '
            WRITE (6,*) ' lt,lr,lmx,length: ',LT,LR,LMX,LENGTH
            MEMALL = (-1)
            RETURN

         END IF

         IF (IPN.EQ.0) THEN
            LAV = LMX - IP - 1

         ELSE IF (IPN.LT.0) THEN
            LAV = -IPN - IP - 2
         END IF

         IF (LAV.GT.LR) THEN
            MEMRY(IP) = LT + 1
            MEMRY(LT+1) = IPN
            IF (IPN.EQ.0) THEN
               MEMRY(2) = LT + 1
            END IF

            IB = IP + 2
            DO 2670 J = 1,LR
               MEMRY(IP+1+J) = 0
 2670       CONTINUE

         ELSE IF (LAV.EQ.LR) THEN
            MEMRY(IP) = -IPN
            IB = IP + 2
            DO 2690 J = 1,LR
               MEMRY(IP+1+J) = 0
 2690       CONTINUE
         END IF

      END IF

      IP = ABS(IPN)
      GO TO 2650

 2660 CONTINUE
      MEMALL = (IB)
      RETURN

      END
