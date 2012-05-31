C $Id: close4.f,v 1.1 1990/11/30 11:04:50 gdsjaar Exp $
C $Log: close4.f,v $
C Revision 1.1  1990/11/30 11:04:50  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]CLOSE4.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CLOSE4 (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   N0, N1, N2, N3, KKK, ERR)
C***********************************************************************
C
C  SUBROUTINE CLOSE4 = CLOSES THE AREA AROUND A FOUR NODE ELEMENT
C
C***********************************************************************
C
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND), LNODES (MLN, MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
C
      LOGICAL ERR
C
C  SET ALL THE LOOP NODES TO BE INTERIOR
C
      LNODES (4, N0) = - 2
      LNODES (4, N1) = - 2
      LNODES (4, N2) = - 2
      LNODES (4, N3) = - 2
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   N0, ERR)
      IF (ERR) GOTO 100
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   N1, ERR)
      IF (ERR) GOTO 100
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   N2, ERR)
      IF (ERR) GOTO 100
      CALL MARKSM (MXND, MLN, LXK, KXL, NXL, LXN, LNODES,
     &   N3, ERR)
      IF (ERR) GOTO 100
C
C  MAKE THE LXK AND THE KXL ARRAY
C
      KKK = KKK+1
      LXK (1, KKK) = LNODES (5, N0)
      LXK (2, KKK) = LNODES (5, N1)
      LXK (3, KKK) = LNODES (5, N2)
      LXK (4, KKK) = LNODES (5, N3)
C
      CALL ADDKXL (MXND, KXL, KKK, LNODES (5, N0))
      CALL ADDKXL (MXND, KXL, KKK, LNODES (5, N1))
      CALL ADDKXL (MXND, KXL, KKK, LNODES (5, N2))
      CALL ADDKXL (MXND, KXL, KKK, LNODES (5, N3))
C
  100 CONTINUE
      RETURN
C
      END
