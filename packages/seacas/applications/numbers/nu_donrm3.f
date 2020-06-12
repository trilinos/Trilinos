C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    See packages/seacas/LICENSE for details

C $Id: donrm3.f,v 1.1 1991/02/21 15:42:57 gdsjaar Exp $
C $Log: donrm3.f,v $
C Revision 1.1  1991/02/21 15:42:57  gdsjaar
C Initial revision
C
      SUBROUTINE DONRM3 (COORD, LTNESS, MAP, DIRCOS, TEMP,
     *    NSEG, NUMNIQ, NUMNP)
      DIMENSION COORD(NUMNP, *), LTNESS(4,*), MAP(*), DIRCOS(5,*),
     *    TEMP(3,*)
C
      DO 10 ISEG = 1, NSEG
          XI = COORD(LTNESS(1,ISEG),1 )
          YI = COORD(LTNESS(1,ISEG),2 )
          ZI = COORD(LTNESS(1,ISEG),3 )
C
          XJ = COORD(LTNESS(2,ISEG),1 )
          YJ = COORD(LTNESS(2,ISEG),2 )
          ZJ = COORD(LTNESS(2,ISEG),3 )
C
          XK = COORD(LTNESS(3,ISEG),1 )
          YK = COORD(LTNESS(3,ISEG),2 )
          ZK = COORD(LTNESS(3,ISEG),3 )
C
          XL = COORD(LTNESS(4,ISEG),1 )
          YL = COORD(LTNESS(4,ISEG),2 )
          ZL = COORD(LTNESS(4,ISEG),3 )
C
          AI =  (YK - YI) * (ZL - ZJ) - (ZK - ZI) * (YL - YJ)
          BJ =  (ZK - ZI) * (XL - XJ) - (XK - XI) * (ZL - ZJ)
          CK =  (XK - XI) * (YL - YJ) - (YK - YI) * (XL - XJ)
          RMAG = SQRT ( AI**2 + BJ**2 + CK**2)
C
          TEMP(1,ISEG) = AI / RMAG
          TEMP(2,ISEG) = BJ / RMAG
          TEMP(3,ISEG) = CK / RMAG
C
   10 CONTINUE
C
      DO 20 I=1,NUMNIQ
          DIRCOS(1,I) = 0.0
          DIRCOS(2,I) = 0.0
          DIRCOS(3,I) = 0.0
   20 CONTINUE
C
      DO 40 ISEG = 1, NSEG
          DO 30 J = 1, 4
              MISEG = MAP( 4 * (ISEG-1) + J )
              DIRCOS(1,MISEG) = DIRCOS(1,MISEG) + TEMP(1,ISEG)
              DIRCOS(2,MISEG) = DIRCOS(2,MISEG) + TEMP(2,ISEG)
              DIRCOS(3,MISEG) = DIRCOS(3,MISEG) + TEMP(3,ISEG)
   30     CONTINUE
   40 CONTINUE
C
C ... NORMALIZE ALL DIRECTION COSINES
C
      DO 50 I = 1, NUMNIQ
          A = DIRCOS(1,I)
          B = DIRCOS(2,I)
          C = DIRCOS(3,I)
          R = SQRT(A**2 + B**2 + C**2)
          DIRCOS(1,I) = A/R
          DIRCOS(2,I) = B/R
          DIRCOS(3,I) = C/R
   50 CONTINUE
      RETURN
      END
