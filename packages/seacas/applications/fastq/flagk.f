C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE FLAGK (NPELEM, NNXK, NXK, MAPDXG, I1, I2, SETFLG, OLD)
C************************************************************************

C  SUBROUTINE FLAGK = FLAGS ELEMENTS FOR PLOTTING OR NOT PLOTTING

C***********************************************************************

C  VARIABLES USED:
C     NPELEMS = NUMBER OF PROCESSED ELEMENTS
C     NXK     = NODES PER ELEMENT ARRAY (CONNECTIVITY)
C     I1      = BEGINNING ELEMENT TO BE FLAGGED
C     I2      = ENDING ELEMENT TO BE FLAGGED
C     SETFLG  = .TRUE. IF THE ELEMENT IS TO BE FLAGGED FOR PLOTTING
C     OLD     = .TRUE. IF THE OLD ELEMENT NUMBERS ARE TO BE USED

C***********************************************************************

C  NOTE:
C     THE ELEMENT IS FLAGGED FOR PLOTTING BY FORCING THE FIRST NODE TO
C     BE POSITIVE AND VICE VERSUS

C***********************************************************************

      DIMENSION NXK(NNXK,NPELEM), MAPDXG(NPELEM)

      LOGICAL SETFLG, OLD

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

      END
