C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE LTNEW (MDIM, LINK)
C***********************************************************************

C  SUBROUTINE LTNEW = LOOKUP TABLE CLEARING FOR DATA POINTER ARRAYS

C***********************************************************************

C  VARIABLES USED:
C     MDIM   = DIMENSION OF LINK ARRAY, AND BASE FOR LOOKUP START
C     LINK   = LOOKUP TABLE ARRAY OF ID'S AND POINTERS
C              LINK(1,I) = ID VALUE STORED IN I'TH ROW (0 IF EMPTY)
C              LINK(2,I) = DATA POINTER ASSOCIATED W/THIS I'TH ID VALUE

C***********************************************************************

      DIMENSION LINK(2,MDIM)

C  ZERO OUT THE ID ROW ONLY

      DO 100 I = 1, MDIM
         LINK(1,I) = 0
  100 CONTINUE

      RETURN

      END
