C $Id: d2grid.f,v 1.1 1990/11/30 11:05:37 gdsjaar Exp $
C $Log: d2grid.f,v $
C Revision 1.1  1990/11/30 11:05:37  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]D2GRID.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE D2GRID (X1, Y1, X2, Y2)
C***********************************************************************
C
C  SUBROUTINE D2GRID = DRAWS A LINE BETWEEN TWO GRIDS
C
C***********************************************************************
C
      DIMENSION X (2), Y (2)
C
      X (1) = X1
      X (2) = X2
      Y (1) = Y1
      Y (2) = Y2
      CALL MPD2VC (1, X (1), Y (1), X (2), Y (2))
      RETURN
C
      END
