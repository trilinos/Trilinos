C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      FUNCTION INDX (N, L, IVAL)
C************************************************************************

C  FUNCTION INDX = FINDS THE INDEX IN L OF IVAL

C************************************************************************

C  NOTE:
C     L MUST BE IN INCREASING ORDER
C     IF IVAL IS NOT IN L,  INDEX=0 IS RETURNED

C***********************************************************************

      DIMENSION L (N)

C  BISECTION SEARCH

      IF (N .LT. 1) THEN
         INDX=0
         RETURN
      ENDIF
      ILO=1
      IHI=N
  100 CONTINUE
      IMID= (ILO + IHI) / 2

C  CONVERGENCE

      IF (IMID .EQ. ILO) THEN
         IF (IVAL .EQ. L (IMID)) THEN
            INDX=IMID
            RETURN
         ELSEIF (IVAL .NE. L (IHI)) THEN
            INDX=0
            RETURN
         ENDIF
         INDX=IHI
         RETURN
      ENDIF

      IF (IVAL .LT. L (IMID)) THEN
         IHI=IMID
      ELSEIF (IVAL .EQ. L (IMID)) THEN
         INDX=IMID
         RETURN
      ELSE
         ILO=IMID
      ENDIF
      GOTO 100

      END
