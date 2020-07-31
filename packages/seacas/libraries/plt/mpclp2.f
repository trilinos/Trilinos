C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION MPCLP2(N,C1X,C1Y,C2X,C2Y)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      CHARACTER*6 SUBNAM
      DIMENSION C1X(*),C1Y(*),C2X(*),C2Y(*)
      PARAMETER (SUBNAM='MPCLP2')

      MPCLP2 = .FALSE.
      IF (N.GT.10) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *               'Too many clipping lines specified; max is 10',2)
         RETURN

      END IF

      IF (N.LT.0) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *               'You cannot specify less than zero clipping lines',
     *               2)
         RETURN

      END IF

      MPCLP2 = .TRUE.
      NCPLIN = N
      DO 2420 I = 1,N
         CPLINE(1,1,I) = C1X(I)
         CPLINE(1,2,I) = C1Y(I)
         CPLINE(2,1,I) = C2X(I)
         CPLINE(2,2,I) = C2Y(I)
 2420 CONTINUE
      RETURN

      END
