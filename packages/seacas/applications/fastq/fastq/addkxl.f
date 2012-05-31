C $Id: addkxl.f,v 1.1 1990/11/30 11:02:54 gdsjaar Exp $
C $Log: addkxl.f,v $
C Revision 1.1  1990/11/30 11:02:54  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]ADDKXL.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE ADDKXL (MXND, KXL, K, L)
C***********************************************************************
C
C  SUBROUTINE ADDKXL = ADDS TO THE LIST OF ELEMENTS FOR THIS LINE
C
C***********************************************************************
C
      DIMENSION KXL (2, 3*MXND)
C
      IF ( KXL(1, L) .EQ. 0) THEN
         KXL(1, L) = K
      ELSE
         KXL(2, L) = K
      ENDIF
      RETURN
C
      END
