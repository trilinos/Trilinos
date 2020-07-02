C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    See packages/seacas/LICENSE for details

C $Id: mesage.f,v 1.1 1990/11/30 11:11:59 gdsjaar Exp $
C $Log: mesage.f,v $
C Revision 1.1  1990/11/30 11:11:59  gdsjaar
C Initial revision
C
C
      SUBROUTINE MESAGE (PROMPT)
C***********************************************************************
C
C  SUBROUTINE MESAGE = PRINTS A MESSAGE ONTO THE SCREEN
C
C***********************************************************************
C
      CHARACTER * (*) PROMPT
C
      IF (PROMPT .EQ. ' ') THEN
         WRITE (*, 10000)
      ELSE
         WRITE (*, 10010)PROMPT
      ENDIF
      RETURN
C
10000 FORMAT ( / )
10010 FORMAT (' ', A)
      END
