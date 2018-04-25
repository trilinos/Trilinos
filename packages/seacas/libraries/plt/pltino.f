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

C $Id: pltino.f,v 1.2 2000/10/25 13:32:35 gdsjaar Exp $ 
C $Log: pltino.f,v $
C Revision 1.2  2000/10/25 13:32:35  gdsjaar
C Modified intrinsic functions to use generic versions to avoid warnings on SGI 64-bit compiles
C
C Revision 1.1  1993/07/16 16:48:25  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      SUBROUTINE PLTINO(MIN,MAX,START,REND,INTER,EXP,NMIN)
      REAL MIN,MAX,INTER,NINTER
      INTEGER EXP
      REAL NI
      DATA SMALL/1.E-4/,DESINT/5./

      DELTA = MAX - MIN
      TMAX = MAX
      TMIN = MIN
      IF (DELTA.NE.0) THEN
         EXP = NINT(LOG10(ABS(DELTA))) - 1

      ELSE
         IF (MIN.NE.0) THEN
            EXP = NINT(LOG10(ABS(MIN))) - 1
            EPS = .01*ABS(MIN)

         ELSE
            EXP = 0
            EPS = .1
         END IF

         TMAX = MAX + EPS
         TMIN = TMIN - EPS
      END IF

      TENEXP = 10.**EXP
      RMIN = 1.E20
      J = 1
      DO 2470 I = 1,5
         NINTER = DELTA/ (FLOAT(I)*TENEXP)
         TEMP = ABS(DESINT-NINTER)
         IF (TEMP.LT.RMIN) THEN
            J = I
            RMIN = TEMP
         END IF

 2470 CONTINUE
      INTER = FLOAT(J)
      IF (DELTA.EQ.0.) THEN
         INTER = 1.
      END IF

      IF (INTER.EQ.1.) THEN
         NMIN = 5

      ELSE IF (INTER.EQ.2.) THEN
         NMIN = 4

      ELSE IF (INTER.EQ.3.) THEN
         NMIN = 6

      ELSE IF (INTER.EQ.4.) THEN
         NMIN = 4

      ELSE IF (INTER.EQ.5.) THEN
         NMIN = 5
      END IF

      TENI = INTER*TENEXP
      IF (TMIN.GE.0.) THEN
         ADD = 0.

      ELSE
         ADD = -1.
      END IF

      RNI = TMIN/TENI
      NI = INT(RNI)
      IF (ABS(RNI-NI).LT.SMALL) THEN
         ADD = 0.
      END IF

      START = (NI+ADD)*TENI
      IF (TMAX.GE.0.) THEN
         ADD = 1.

      ELSE
         ADD = 0.
      END IF

      RNI = TMAX/TENI
      NI = INT(RNI)
      IF (ABS(RNI-NI).LT.SMALL) THEN
         ADD = 0.
      END IF

      REND = (NI+ADD)*TENI
      START = START/TENEXP
      REND = REND/TENEXP
      IF (REND.NE.0.) THEN
 2490    IF (.NOT. (ABS(REND).GT.10.)) GO TO 2500
         REND = REND/10.
         START = START/10.
         INTER = INTER/10.
         EXP = EXP + 1
         GO TO 2490

 2500    CONTINUE
      END IF

      IF (START.NE.0.) THEN
 2510    IF (.NOT. (ABS(START).LT.1.)) GO TO 2520
         REND = REND*10.
         START = START*10.
         INTER = INTER*10.
         EXP = EXP - 1
         GO TO 2510

 2520    CONTINUE
      END IF

      IF (START.EQ.0 .OR. TMIN.EQ.0) THEN
         RETURN

      END IF

      IF (START.EQ.REND) THEN
         REND = START + INTER
      END IF

      IF (ABS(START-REND).EQ.INTER) THEN
         NMIN = 10
      END IF

      RETURN

      END
