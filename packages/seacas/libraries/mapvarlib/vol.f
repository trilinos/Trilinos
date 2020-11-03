C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
*DECK,VOL
      SUBROUTINE VOL (ITYPE,XX,YY,ZZ,VOLUME)

C     ******************************************************************

C     SUBROUTINE TO FIND THE VOLUME OF AN ELEMENT

C     Called by ELTON3, ELTOND

C     ******************************************************************

      DIMENSION XX(*), YY(*), ZZ(*)

C     ******************************************************************

      IF (ITYPE .EQ. 3) THEN

C 4-node quad

        A1 = XX(1) - XX(3)
        B1 = YY(1) - YY(3)
        A2 = XX(2) - XX(4)
        B2 = YY(2) - YY(4)
        CCR = (A1 * B2) - (A2 * B1)
        VOLUME = 0.5 * ABS(CCR)

        RETURN

      ELSE IF (ITYPE .EQ. 10) THEN

C 8-node hex

        GR1 = ( YY(2)*(ZZ(6)-ZZ(3)-ZZ(4)+ZZ(5)) + YY(3)*(ZZ(2)-ZZ(4))
     &        + YY(4)*(ZZ(3)-ZZ(8)-ZZ(5)+ZZ(2))
     &        + YY(5)*(ZZ(8)-ZZ(6)-ZZ(2)+ZZ(4))
     &        + YY(6)*(ZZ(5)-ZZ(2)) + YY(8)*(ZZ(4)-ZZ(5)) ) / 12.

        GR2 = ( YY(3)*(ZZ(7)-ZZ(4)-ZZ(1)+ZZ(6)) + YY(4)*(ZZ(3)-ZZ(1))
     &        + YY(1)*(ZZ(4)-ZZ(5)-ZZ(6)+ZZ(3))
     &        + YY(6)*(ZZ(5)-ZZ(7)-ZZ(3)+ZZ(1))
     &        + YY(7)*(ZZ(6)-ZZ(3)) + YY(5)*(ZZ(1)-ZZ(6)) ) / 12.

        GR3 = ( YY(4)*(ZZ(8)-ZZ(1)-ZZ(2)+ZZ(7)) + YY(1)*(ZZ(4)-ZZ(2))
     &        + YY(2)*(ZZ(1)-ZZ(6)-ZZ(7)+ZZ(4))
     &        + YY(7)*(ZZ(6)-ZZ(8)-ZZ(4)+ZZ(2))
     &        + YY(8)*(ZZ(7)-ZZ(4)) + YY(6)*(ZZ(2)-ZZ(7)) ) / 12.

        GR4 = ( YY(1)*(ZZ(5)-ZZ(2)-ZZ(3)+ZZ(8)) + YY(2)*(ZZ(1)-ZZ(3))
     &        + YY(3)*(ZZ(2)-ZZ(7)-ZZ(8)+ZZ(1))
     &        + YY(8)*(ZZ(7)-ZZ(5)-ZZ(1)+ZZ(3))
     &        + YY(5)*(ZZ(8)-ZZ(1)) + YY(7)*(ZZ(3)-ZZ(8)) ) / 12.

        GR5 = ( YY(8)*(ZZ(4)-ZZ(7)-ZZ(6)+ZZ(1)) + YY(7)*(ZZ(8)-ZZ(6))
     &        + YY(6)*(ZZ(7)-ZZ(2)-ZZ(1)+ZZ(8))
     &        + YY(1)*(ZZ(2)-ZZ(4)-ZZ(8)+ZZ(6))
     &        + YY(4)*(ZZ(1)-ZZ(8)) + YY(2)*(ZZ(6)-ZZ(1)) ) / 12.

        GR6 = ( YY(5)*(ZZ(1)-ZZ(8)-ZZ(7)+ZZ(2)) + YY(8)*(ZZ(5)-ZZ(7))
     &        + YY(7)*(ZZ(8)-ZZ(3)-ZZ(2)+ZZ(5))
     &        + YY(2)*(ZZ(3)-ZZ(1)-ZZ(5)+ZZ(7))
     &        + YY(1)*(ZZ(2)-ZZ(5)) + YY(3)*(ZZ(7)-ZZ(2)) ) / 12.

        GR7 = ( YY(6)*(ZZ(2)-ZZ(5)-ZZ(8)+ZZ(3)) + YY(5)*(ZZ(6)-ZZ(8))
     &        + YY(8)*(ZZ(5)-ZZ(4)-ZZ(3)+ZZ(6))
     &        + YY(3)*(ZZ(4)-ZZ(2)-ZZ(6)+ZZ(8))
     &        + YY(2)*(ZZ(3)-ZZ(6)) + YY(4)*(ZZ(8)-ZZ(3)) ) / 12.

        GR8 = ( YY(7)*(ZZ(3)-ZZ(6)-ZZ(5)+ZZ(4)) + YY(6)*(ZZ(7)-ZZ(5))
     &        + YY(5)*(ZZ(6)-ZZ(1)-ZZ(4)+ZZ(7))
     &        + YY(4)*(ZZ(1)-ZZ(3)-ZZ(7)+ZZ(5))
     &        + YY(3)*(ZZ(4)-ZZ(7)) + YY(1)*(ZZ(5)-ZZ(4)) ) / 12.

        VOLUME = XX(1) * GR1 + XX(2) * GR2 + XX(3) * GR3 + XX(4) * GR4
     &         + XX(5) * GR5 + XX(6) * GR6 + XX(7) * GR7 + XX(8) * GR8

        RETURN

      ELSE IF (ITYPE .EQ. 13) THEN

C 4-node shell, ignore curvature

        A1 = XX(1) - XX(3)
        B1 = YY(1) - YY(3)
        C1 = ZZ(1) - ZZ(3)
        A2 = XX(2) - XX(4)
        B2 = YY(2) - YY(4)
        C2 = ZZ(2) - ZZ(4)

        ACR = (B1 * C2) - (B2 * C1)
        BCR = (A1 * C2) - (A2 * C1)
        CCR = (A1 * B2) - (A2 * B1)

        VOLUME = 0.5 * SQRT(ACR*ACR + BCR*BCR + CCR*CCR)

        RETURN
      ELSE IF (ITYPE .EQ. 6) THEN

C 4-node tet
C (process Key's 8-node tet the same way using first 4 nodes)

        Y12 = YY(1) - YY(2)
        Y13 = YY(1) - YY(3)
        Y14 = YY(1) - YY(4)
        Y24 = YY(2) - YY(4)
        Y34 = YY(3) - YY(4)
        Z12 = ZZ(1) - ZZ(2)
        Z13 = ZZ(1) - ZZ(3)
        Z14 = ZZ(1) - ZZ(4)
        Z24 = ZZ(2) - ZZ(4)
        Z34 = ZZ(3) - ZZ(4)

        BX1 = (Y34*Z24 - Y24*Z34) / 6.0
        BX2 = (Y13*Z14 - Y14*Z13) / 6.0
        BX3 = (Y14*Z12 - Y12*Z14) / 6.0
        BX4 = (Y12*Z13 - Y13*Z12) / 6.0

        VOLUME = BX1*XX(1) + BX2*XX(2) + BX3*XX(3) + BX4*XX(4)
      ELSE
        CALL ERROR('VOL','ELEMENT TYPE',' ',ITYPE,
     &             'NOT YET IMPLEMENTED',0,' ',' ',1)
      END IF
      RETURN
      END
