C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    See packages/seacas/LICENSE for details

C $Id: erasec.f,v 1.1 1990/11/30 11:06:45 gdsjaar Exp $
C $Log: erasec.f,v $
C Revision 1.1  1990/11/30 11:06:45  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]ERASEC.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE ERASEC (OLDCUR)
C***********************************************************************
C
C  SUBROUTINE ERASEC = DEACTIVATES THE CROSSHAIRS
C
C***********************************************************************
C
      LOGICAL OLDCUR
C
      IF (OLDCUR) THEN
         WRITE (*,*) CHAR(27)//'G0'
         OLDCUR = .FALSE.
      ENDIF
      RETURN
C
      END
