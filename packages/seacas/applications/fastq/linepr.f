C $Id: linepr.f,v 1.1 1990/11/30 11:11:02 gdsjaar Exp $
C $Log: linepr.f,v $
C Revision 1.1  1990/11/30 11:11:02  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]LINEPR.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE LINEPR (ML, MP, LINKP, LCON, II, I1, I2, I3, J1, J2,
     &   J3)
C***********************************************************************
C
C  SUBROUTINE LINEPR = GETS THE LINE PARAMETERS
C***********************************************************************
C
      DIMENSION LCON(3, ML)
      DIMENSION LINKP(2, MP)
C
      LOGICAL ADDLNK
C
      ADDLNK = .FALSE.
      I1 = LCON (1, II)
      I2 = LCON (2, II)
      I3 = LCON (3, II)
      CALL LTSORT (MP, LINKP, I1, J1, ADDLNK)
      CALL LTSORT (MP, LINKP, I2, J2, ADDLNK)
      IF (I3 .NE. 0) THEN
         CALL LTSORT (MP, LINKP, IABS(I3), J3, ADDLNK)
      ELSE
         J3 = 0
      ENDIF
C
      RETURN
C
      END
