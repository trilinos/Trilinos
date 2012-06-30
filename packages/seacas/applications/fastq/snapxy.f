C $Id: snapxy.f,v 1.1 1990/11/30 11:15:59 gdsjaar Exp $
C $Log: snapxy.f,v $
C Revision 1.1  1990/11/30 11:15:59  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]SNAPXY.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE SNAPXY (MP, MSNAP, N, IPOINT, COOR, LINKP, SNAPDX,
     &   NSNAP)
C***********************************************************************
C
C  SUBROUTINE SNAPXY = ADDS SNAP GRIDS AT EVERY X, Y LOCATION
C
C***********************************************************************
C
      DIMENSION IPOINT (MP), COOR (2, MP), LINKP (2, MP)
      DIMENSION SNAPDX (2, MSNAP), NSNAP (2)
C
      LOGICAL ERR, ADDLNK
C
      ADDLNK = .FALSE.
C
      DO 100 I = 1, N
         CALL LTSORT (MP, LINKP, IABS (IPOINT (I)), II, ADDLNK)
         IF (II.GT.0)THEN
            INDEX = 1
            CALL ADDSNP (MSNAP, SNAPDX, NSNAP, INDEX, COOR (1, II), ERR)
            IF  (ERR) RETURN
            INDEX = 2
            CALL ADDSNP (MSNAP, SNAPDX, NSNAP, INDEX, COOR (2, II), ERR)
            IF  (ERR) RETURN
         ENDIF
  100 CONTINUE
C
      RETURN
C
      END
