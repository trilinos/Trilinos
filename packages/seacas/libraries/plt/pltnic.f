C Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
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

C $Id: pltnic.f,v 1.2 2000/10/25 13:32:35 gdsjaar Exp $ 
C $Log: pltnic.f,v $
C Revision 1.2  2000/10/25 13:32:35  gdsjaar
C Modified intrinsic functions to use generic versions to avoid warnings on SGI 64-bit compiles
C
C Revision 1.1  1993/07/16 16:49:01  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      SUBROUTINE PLTNIC(X,TYPE,FN,NE,INTER,NMIN)
      CHARACTER*(*) TYPE
      CHARACTER*1 TTYPE
      REAL FNICE(11),INTERA(11),INTER
      INTEGER NMINA(11)
      DATA FNICE/0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10./
      DATA INTERA/1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1./
      DATA NMINA/10,10,10,10,5,5,5,5,5,5,5/

      TTYPE = TYPE
      CALL CHRUP(TTYPE,TTYPE)
      IF (X.EQ.0.) THEN
         FN = 0.
         NE = 0
         RETURN

      END IF

      XTEMP = ABS(X)
      IF (X.LT.0. .AND. TTYPE.EQ.'O') THEN
         TTYPE = 'U'

      ELSE IF (X.LT.0. .AND. TTYPE.EQ.'U') THEN
         TTYPE = 'O'
      END IF

      E1 = LOG10(XTEMP)
      JE = INT(E1)
      IF (JE.LE.0. .AND. XTEMP.LT.1.) THEN
         JE = JE - 1
      END IF

      IF (JE.LT.0) THEN
         E2 = 1./10.** (-JE)

      ELSE IF (JE.EQ.0) THEN
         E2 = 1.

      ELSE IF (JE.GT.0) THEN
         E2 = 10.**JE
      END IF

      F1 = XTEMP/E2
      IF (TTYPE.EQ.'O') THEN
         DO 2910 I = 1,11
            IF (F1/1.007.LE.FNICE(I)) THEN
               GO TO 2920

            END IF

 2910    CONTINUE
 2920    CONTINUE
      END IF

      IF (TTYPE.EQ.'U') THEN
         DO 2930 I = 11,1,-1
            IF (F1*1.007.GE.FNICE(I)) THEN
               GO TO 2940

            END IF

 2930    CONTINUE
 2940    CONTINUE
      END IF

      IF (I.GT.11) THEN
         I = 11
      END IF

      IF (I.LT.1) THEN
         I = 1
      END IF

      FN = FNICE(I)
      INTER = INTERA(I)
      NMIN = NMINA(I)
      NE = JE
      IF (X.LT.0.) THEN
         FN = -FN
      END IF

      IF (FN.EQ.0.) THEN
         NE = 0
      END IF

      RETURN

      END
