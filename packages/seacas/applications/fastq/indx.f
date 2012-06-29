C $Id: indx.f,v 1.1 1990/11/30 11:09:32 gdsjaar Exp $
C $Log: indx.f,v $
C Revision 1.1  1990/11/30 11:09:32  gdsjaar
C Initial revision
C
C
CC* FILE: [.RENUM]INDX.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      FUNCTION INDX (N, L, IVAL)
C************************************************************************
C
C  FUNCTION INDX = FINDS THE INDEX IN L OF IVAL
C
C************************************************************************
C
C  NOTE:
C     L MUST BE IN INCREASING ORDER
C     IF IVAL IS NOT IN L,  INDEX=0 IS RETURNED
C
C***********************************************************************
C
      DIMENSION L (N)
C
C  BISECTION SEARCH
C
      IF (N .LT. 1) THEN
         INDX=0
         RETURN
      ENDIF
      ILO=1
      IHI=N
  100 CONTINUE
      IMID= (ILO + IHI) / 2
C
C  CONVERGENCE
C
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
C
      IF (IVAL .LT. L (IMID)) THEN
         IHI=IMID
      ELSEIF (IVAL .EQ. L (IMID)) THEN
         INDX=IMID
         RETURN
      ELSE
         ILO=IMID
      ENDIF
      GOTO 100
C
      END
