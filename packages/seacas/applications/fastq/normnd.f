C $Id: normnd.f,v 1.1 1990/11/30 11:12:49 gdsjaar Exp $
C $Log: normnd.f,v $
C Revision 1.1  1990/11/30 11:12:49  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]NORMND.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE NORMND (NPNODE, BMESUR, RMAX)
C***********************************************************************
C
C  SUBROUTINE NORMND = NORMALIZES A NODE VARIABLE
C
C***********************************************************************
C
      DIMENSION BMESUR(NPNODE)
C
      BMIN = BMESUR(1)
      BMAX = BMESUR(1)
      DO 100 NODE = 2, NPNODE
         BMAX = AMAX1 (BMESUR(NODE), BMAX)
         BMIN = AMIN1 (BMESUR(NODE), BMIN)
  100 CONTINUE
C
      BMAX = BMAX - BMIN
      DO 110 NODE = 1, NPNODE
         BMESUR(NODE) = BMESUR(NODE) - BMIN
  110 CONTINUE
C
C  RMAX = MAXIMUM RATIO FOR PLATEAU VALUES
C
      DO 120 NODE = 1, NPNODE
         IF (BMESUR (NODE) .GE. (BMAX * RMAX)) THEN
            BMESUR(NODE) = 1.0
         ELSE
            BMESUR (NODE) = BMESUR(NODE) / (BMAX * RMAX)
         ENDIF
  120 CONTINUE
C
      RETURN
C
      END
