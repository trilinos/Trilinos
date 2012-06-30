C $Id: jumplp.f,v 1.1 1990/11/30 11:10:38 gdsjaar Exp $
C $Log: jumplp.f,v $
C Revision 1.1  1990/11/30 11:10:38  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]JUMPLP.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      FUNCTION JUMPLP (MXND, MLN, LNODES, INOW, IJUMP)
C***********************************************************************
C
C  FUNCTION JUMPLP = JUMPS IJUMP STEPS FORWARD AROUND THE CLOSED LOOP
C
C***********************************************************************
C
      DIMENSION LNODES (MLN, MXND)
C
      JUMPLP = INOW
      DO 100 I = 1, IJUMP
         JUMPLP = LNODES (3, JUMPLP)
  100 CONTINUE
      RETURN
C
      END
