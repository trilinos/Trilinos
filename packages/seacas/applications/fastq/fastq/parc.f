C $Id: parc.f,v 1.1 1990/11/30 11:13:03 gdsjaar Exp $
C $Log: parc.f,v $
C Revision 1.1  1990/11/30 11:13:03  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]PARC.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      FUNCTION PARC (AL, TCOEF)
C***********************************************************************
C
C    SUBROUTINE PARC = CALCULATES PARABOLIC ARC LOCATIONS
C
C***********************************************************************
C
      PARC = 0.5 * (SQRT (1.0 + AL **2) * AL +
     &   ALOG (SQRT (1.0 + AL **2) + AL)) / TCOEF
C
      END
