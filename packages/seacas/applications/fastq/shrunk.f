C $Id: shrunk.f,v 1.1 1990/11/30 11:15:38 gdsjaar Exp $
C $Log: shrunk.f,v $
C Revision 1.1  1990/11/30 11:15:38  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]SHRUNK.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      LOGICAL FUNCTION SHRUNK (RATIO, NROW)
C***********************************************************************
C
C  FUNCTION SHRUNK = LOGICAL FUNCTION THAT RETURNS TRUE IF THE ELEMENT
C                    SIZE IS DIMINISHING WITH ROW DEPTH
C
C***********************************************************************
C
      DATA TOLER1 /.85/, TOLER2 /.75/, TOLER3 /.6/
C
      IF ((NROW .GE. 3) .AND. (RATIO .LT. TOLER1)) THEN
         SHRUNK = .TRUE.
      ELSEIF ((NROW .GE. 2) .AND. (RATIO .LT. TOLER2)) THEN
         SHRUNK = .TRUE.
      ELSEIF ((NROW .GE. 1) .AND. (RATIO .LT. TOLER3)) THEN
         SHRUNK = .TRUE.
      ELSE
         SHRUNK = .FALSE.
      ENDIF
C
      RETURN
      END
