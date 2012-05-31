C $Id: ltnew.f,v 1.1 1990/11/30 11:11:40 gdsjaar Exp $
C $Log: ltnew.f,v $
C Revision 1.1  1990/11/30 11:11:40  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]LTNEW.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE LTNEW (MDIM, LINK)
C***********************************************************************
C
C  SUBROUTINE LTNEW = LOOKUP TABLE CLEARING FOR DATA POINTER ARRAYS
C
C***********************************************************************
C
C  VARIABLES USED:
C     MDIM   = DIMENSION OF LINK ARRAY, AND BASE FOR LOOKUP START
C     LINK   = LOOKUP TABLE ARRAY OF ID'S AND POINTERS
C              LINK(1,I) = ID VALUE STORED IN I'TH ROW (0 IF EMPTY)
C              LINK(2,I) = DATA POINTER ASSOCIATED W/THIS I'TH ID VALUE
C
C***********************************************************************
C
      DIMENSION LINK(2,MDIM)
C
C  ZERO OUT THE ID ROW ONLY
C
      DO 100 I = 1, MDIM
         LINK(1,I) = 0
  100 CONTINUE
C
      RETURN
C
      END
