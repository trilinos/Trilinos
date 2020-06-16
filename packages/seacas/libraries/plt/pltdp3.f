C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C $Id: pltdp3.f,v 1.4 2000/10/25 18:55:02 gdsjaar Exp $
C $Log: pltdp3.f,v $
C Revision 1.4  2000/10/25 18:55:02  gdsjaar
C In the pltli? functions, check for N==0 before doing any array
C accesses.
C
C Also changed all references to 'mask' to be arrays where they were
C scalars since downstream code seems to treat them as arrays.
C
C Revision 1.3  1993/07/19 17:06:38  gdsjaar
C Changed hex constants back to preceding X, --needed on cray. Works
C either way on other systems.
C
c Revision 1.2  1993/07/16  17:33:15  gdsjaar
c Integer constant too big on sun, replaced it with hexadecimal notation
c
c Revision 1.1  1993/07/16  16:48:02  gdsjaar
c Changed plt to library rather than single source file.
c
C=======================================================================
      SUBROUTINE PLTDP3(MAP,N,PX,PY,PZ)
      REAL MAP(*),PX(*),PY(*),PZ(*)
      DIMENSION Q1(3),V1(3),Q2(3),V2(3),TPX(32),TPY(32),QX(32),QY(32)
      INTEGER MASK(1)
      include 'izbit.inc'

      DO 2410 L = 1,3
         V1(L) = MAP(18+L-1) + MAP(15)*MAP(27+L-1)
         Q1(L) = MAP(L+27-1)
         V2(L) = MAP(18+L-1) + MAP(16)*MAP(27+L-1)
         Q2(L) = -MAP(L+27-1)
 2410 CONTINUE
      J = 0
      KM = 0
 2430 IF (.NOT. (J.LT.N)) GO TO 2440
      JN = MIN(N-J,32)
      J1 = J + 1
      J = J + JN
      KM = KM + 1
      MASK(1) = -1
      CALL PLTCP3(JN,MASK,PX(J1),PY(J1),PZ(J1),V1,Q1)
      CALL PLTCP3(JN,MASK,PX(J1),PY(J1),PZ(J1),V2,Q2)
      IF (MAP(17).EQ.1.) THEN
         DO 2450 K = 1,JN
            JB = IZBIT(K)
            IF (IAND(JB,MASK(1)).NE.0) THEN
               PMS = (PX(K+J1-1)-MAP(18))*MAP(27) +
     *               (PY(K+J1-1)-MAP(19))*MAP(28) +
     *               (PZ(K+J1-1)-MAP(20))*MAP(29)
               R = MAP(30)/PMS
               TPX(K) = R* ((PX(K+J1-1)-MAP(18))*MAP(21)+
     *                  (PY(K+J1-1)-MAP(19))*MAP(22)+
     *                  (PZ(K+J1-1)-MAP(20))*MAP(23))
               TPY(K) = R* ((PX(K+J1-1)-MAP(18))*MAP(24)+
     *                  (PY(K+J1-1)-MAP(19))*MAP(25)+
     *                  (PZ(K+J1-1)-MAP(20))*MAP(26))
            END IF

 2450    CONTINUE
      END IF

      IF (MAP(17).EQ.-1.) THEN
      END IF

      CALL PLTMP2(MAP,JN,MASK,TPX,TPY,QX,QY)
      CALL PLTPTM(JN,MASK,QX,QY)
      GO TO 2430

 2440 CONTINUE
      RETURN

      END
