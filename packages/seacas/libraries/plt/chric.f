C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CHRIC(JIN,DS,ND)
      CHARACTER*(*) DS
      CHARACTER SDG*16

      DS = ' '
      J = ABS(JIN)
      ND = 0
 2030 CONTINUE
      ND = ND + 1
      KD = MOD(J,10)
      J = J/10
      SDG(ND:ND) = CHAR(ICHAR('0')+KD)

      IF (.NOT. (J.EQ.0)) GO TO 2030
      IJ = 1
      IF (JIN.LT.0) THEN
         IJ = 2
         DS(1:1) = '-'
         ND = ND + 1
      END IF

      I = IJ

 2060 IF (.NOT. (I.LE.ND)) GO TO 2080
      DS(I:I) = SDG(ND-I+1:ND-I+1)
      I = I + 1
      GO TO 2060

 2080 CONTINUE
      RETURN

      END
