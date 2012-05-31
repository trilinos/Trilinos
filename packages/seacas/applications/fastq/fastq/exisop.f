C $Id: exisop.f,v 1.2 1998/07/14 18:18:48 gdsjaar Exp $
C $Log: exisop.f,v $
C Revision 1.2  1998/07/14 18:18:48  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:07:04  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:07:03  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]EXISOP.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE EXISOP (MXND, XN, YN, LNODES, ANGLE, N1, XNEW, YNEW)
C***********************************************************************
C
C  SUBROUTINE EXISOP = CALCULATES A POSITION TO MAKE A PARALLELPIPED
C
C***********************************************************************
C
      DIMENSION XN (MXND), YN (MXND), LNODES (7, MXND), ANGLE (MXND)
C
      N0 = LNODES (2, N1)
      N2 = LNODES (3, N1)
C
      XNEW = XN (N0) + XN (N2) - XN (N1)
      YNEW = YN (N0) + YN (N2) - YN (N1)
C
      RETURN
C
      END
