C $Id: flagk.f,v 1.1 1990/11/30 11:07:35 gdsjaar Exp $
C $Log: flagk.f,v $
C Revision 1.1  1990/11/30 11:07:35  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]FLAGK.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE FLAGK (NPELEM, NNXK, NXK, MAPDXG, I1, I2, SETFLG, OLD)
C************************************************************************
C
C  SUBROUTINE FLAGK = FLAGS ELEMENTS FOR PLOTTING OR NOT PLOTTING
C
C***********************************************************************
C
C  VARIABLES USED:
C     NPELEMS = NUMBER OF PROCESSED ELEMENTS
C     NXK     = NODES PER ELEMENT ARRAY (CONNECTIVITY)
C     I1      = BEGINNING ELEMENT TO BE FLAGGED
C     I2      = ENDING ELEMENT TO BE FLAGGED
C     SETFLG  = .TRUE. IF THE ELEMENT IS TO BE FLAGGED FOR PLOTTING
C     OLD     = .TRUE. IF THE OLD ELEMENT NUMBERS ARE TO BE USED
C
C***********************************************************************
C
C  NOTE:
C     THE ELEMENT IS FLAGGED FOR PLOTTING BY FORCING THE FIRST NODE TO
C     BE POSITIVE AND VICE VERSUS
C
C***********************************************************************
C
      DIMENSION NXK(NNXK,NPELEM), MAPDXG(NPELEM)
C
      LOGICAL SETFLG, OLD
C
      IF (OLD) THEN
         IF (SETFLG) THEN
            DO 100 I = I1, I2
               NXK(1,I) = IABS (NXK (1,I))
  100       CONTINUE
         ELSE
            DO 110 I = I1, I2
               NXK(1,I) = -IABS (NXK(1,I))
  110       CONTINUE
         ENDIF
      ELSE
         IF (SETFLG) THEN
            DO 120 I = I1, I2
               J = MAPDXG (I)
               NXK(1,J) = IABS (NXK (1,J))
  120       CONTINUE
         ELSE
            DO 130 I = I1, I2
               J = MAPDXG(I)
               NXK(1,J) = -IABS (NXK (1,J))
  130       CONTINUE
         ENDIF
      ENDIF
      RETURN
C
      END
