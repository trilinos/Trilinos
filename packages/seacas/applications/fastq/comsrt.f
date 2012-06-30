C $Id: comsrt.f,v 1.1 1990/11/30 11:05:12 gdsjaar Exp $
C $Log: comsrt.f,v $
C Revision 1.1  1990/11/30 11:05:12  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]COMSRT.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE COMSRT (MXND, MXCORN, MXPICK, MLN, LNODES, LCORN,
     &   NCORN, ICOMB, ITYPE, NPICK)
C***********************************************************************
C
C  SUBROUTINE COMSRT = THIS SUBROUTINE GETS ALL THE COMBINATIONS
C                      POSSIBLE OF CORNERS, SIDES, AND DISSECION
C                      NODES
C
C***********************************************************************
C
C  VARIABLES USED:
C     ICOMB  = THE DIFFERENT VALID PRIMITIVE COMBINATIONS
C                 ICOMB (I, J) WHERE I IS THE COMBINATION NUMBER
C                              AND J IS THE CORNER NUMBER.
C                 THE VALUE OF ICOMB (I,J) IS 1 FOR CORNER AND 0 FOR
C                 FOR A SIDE INTERPRETATION
C     ITYPE  = THE TYPE OF PRIMITIVE OR NUMBER OF CORNERS IN THIS
C              COMBINATION.
C                  = 0 FOR A CIRCLE
C                  = 1 FOR A TEARDROP
C                  = 2 FOR A SEMICIRCLE
C                  = 3 FOR A TRIANGLE
C                  = 4 FOR A RECTANGLE
C                  = >4 OTHERWISE
C
C
C***********************************************************************
C
      DIMENSION LNODES (MLN, MXND), ICOMB (MXCORN, MXPICK)
      DIMENSION ITYPE (MXPICK), LCORN (MXCORN)
C
      NPICK = 1
C
      DO 100 I = 1, MXPICK
         ITYPE (I) = 0
  100 CONTINUE
C
      DO 150 I = 1, NCORN
C
C  PUT PURE CORNER AND CORNER/SIDE DESIGNATIONS IN ORDER
C
         ITEST = LNODES (6, LCORN (I))
         IF (ITEST .LE. 2) THEN
            DO 110 J = 1, NPICK
               ICOMB (I, J) = 1
               ITYPE (J) = ITYPE (J) + 1
  110       CONTINUE
            IF (ITEST .EQ. 2) THEN
               DO 130 J = NPICK + 1, NPICK * 2
                  DO 120 K = 1, I - 1
                     ICOMB (K, J) = ICOMB (K, J - NPICK)
  120             CONTINUE
                  ITYPE (J) = ITYPE (J - NPICK) - 1
                  ICOMB (I, J) = 0
  130          CONTINUE
               NPICK = NPICK * 2
            ENDIF
C
C  PUT PURE SIDE AND SIDE/DISSECTIONS DESIGNATIONS IN ORDER
C
         ELSEIF (ITEST .LE. 4) THEN
            DO 140 J = 1, NPICK
               ICOMB (I, J) = 0
  140       CONTINUE
         ENDIF
C
  150 CONTINUE
C
      RETURN
C
      END
