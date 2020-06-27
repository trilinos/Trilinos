C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: mpmul4.f,v 1.4 1993/07/19 17:06:30 gdsjaar Exp $
C $Log: mpmul4.f,v $
C Revision 1.4  1993/07/19 17:06:30  gdsjaar
C Changed hex constants back to preceding X, --needed on cray. Works
C either way on other systems.
C
c Revision 1.3  1993/07/16  22:11:18  gdsjaar
c Unrolled do loops to speed up execution.
c
c Revision 1.2  1993/07/16  17:33:06  gdsjaar
c Integer constant too big on sun, replaced it with hexadecimal notation
c
c Revision 1.1  1993/07/16  16:47:20  gdsjaar
c Changed plt to library rather than single source file.
c
C=======================================================================
      SUBROUTINE MPMUL4(N,MASK,ARR1,ARR2,ARR3,ARR4,MAT,RES1,RES2,RES3,
     *                  RES4)
      DIMENSION ARR1(*),ARR2(*),ARR3(*),ARR4(*),RES1(*),
     *          RES2(*),RES3(*),RES4(*),MAT(4,4)
      REAL MAT
      include 'izbit.inc'

      IF (MASK.EQ.0) THEN
         RETURN

      END IF

      DO 3160 I = 1,N
         IF (IAND(MASK,IZBIT(I)).NE.0) THEN
           RES1(I) = MAT(1,1)*ARR1(I) + MAT(2,1)*ARR2(I) +
     *       MAT(3,1)*ARR3(I) + MAT(4,1)*ARR4(I)
           RES2(I) = MAT(1,2)*ARR1(I) + MAT(2,2)*ARR2(I) +
     *       MAT(3,2)*ARR3(I) + MAT(4,2)*ARR4(I)
           RES3(I) = MAT(1,3)*ARR1(I) + MAT(2,3)*ARR2(I) +
     *       MAT(3,3)*ARR3(I) + MAT(4,3)*ARR4(I)
           RES4(I) = MAT(1,4)*ARR1(I) + MAT(2,4)*ARR2(I) +
     *       MAT(3,4)*ARR3(I) + MAT(4,4)*ARR4(I)
         END IF

 3160 CONTINUE
      RETURN

      END
