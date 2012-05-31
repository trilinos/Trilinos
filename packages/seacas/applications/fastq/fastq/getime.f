C $Id: getime.f,v 1.1 1990/11/30 11:08:13 gdsjaar Exp $
C $Log: getime.f,v $
C Revision 1.1  1990/11/30 11:08:13  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]FRSTRM.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GETIME (TIME)
C***********************************************************************
C
C  SUBROUTINE GETIME = GETS THE CPU TIME USED BY THE CURRENT PROCESS
C
C***********************************************************************
C
      CALL EXCPUS (TIME)
      RETURN
      END
