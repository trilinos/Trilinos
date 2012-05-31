C $Id: d2node.f,v 1.1 1990/11/30 11:05:40 gdsjaar Exp $
C $Log: d2node.f,v $
C Revision 1.1  1990/11/30 11:05:40  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]D2NODE.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE D2NODE (MXND, XN, YN, NODE1, NODE2)
C***********************************************************************
C
C  SUBROUTINE D2NODE = DRAWS A LINE BETWEEN TWO NODES
C
C***********************************************************************
C
      DIMENSION X (2), Y (2), XN (MXND), YN (MXND)
C
      X (1) = XN (NODE1)
      X (2) = XN (NODE2)
      Y (1) = YN (NODE1)
      Y (2) = YN (NODE2)
      CALL MPD2VC (1, X (1), Y (1), X (2), Y (2))
      RETURN
C
      END
