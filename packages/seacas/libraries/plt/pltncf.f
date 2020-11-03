C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PLTNCF(X,TYPE,FN,NE)
      CHARACTER*(*) TYPE
      CHARACTER*1 TTYPE
      REAL FNICE(17)
      DATA FNICE/-10.,-8.,-6.,-5.,-4.,-3.,-2.,-1.,0.,1.,2.,3.,4.,5.,6.,
     *     8.,10./

      TTYPE = TYPE
      CALL CHRUP(TTYPE,TTYPE)
      IF (X.EQ.0.) THEN
         FN = 0.
         RETURN

      END IF

      F1 = X/10.**NE
      IF (TTYPE.EQ.'O') THEN
         DO 2870 I = 1,17
            IF (F1.LE.FNICE(I)) THEN
               GO TO 2880

            END IF

 2870    CONTINUE
 2880    CONTINUE
      END IF

      IF (TTYPE.EQ.'U') THEN
         DO 2890 I = 17,1,-1
            IF (F1.GE.FNICE(I)) THEN
               GO TO 2900

            END IF

 2890    CONTINUE
 2900    CONTINUE
      END IF

      IF (I.GT.17) THEN
         I = 17
      END IF

      IF (I.LT.1) THEN
         I = 1
      END IF

      FN = FNICE(I)
      RETURN

      END
