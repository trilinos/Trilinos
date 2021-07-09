C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

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
