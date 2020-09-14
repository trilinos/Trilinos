C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE SETAXS (X, Y)
C***********************************************************************

C  SETAXS = SETS UP THE AXIS AS NEEDED FOR PLOTTING

C***********************************************************************

      DIMENSION BUFF (11), X (2), Y (2)

C  GET THE AXIS ON THE CORRECT LOCATION ON THE GRAPH

      CALL MP2PT (1, X (1), Y (1), X01, Y01, MASK)
      CALL MP2PT (1, X (2), Y (2), X02, Y02, MASK)

C  FORCE X AND Y LIMITS ON THE GRAPH

      BUFF(1) = X01
      CALL PLTSTG (1, BUFF)
      BUFF(1) = Y01
      CALL PLTSTG (2, BUFF)
      BUFF(1) = X02 - X01
      CALL PLTSTG (3, BUFF)
      BUFF(1) = Y02 - Y01
      CALL PLTSTG (4, BUFF)

C  TURN OFF THE ZERO LINE PLOT

      BUFF(1) = 0.
      CALL PLTSTG (37, BUFF)

C  GET NICE INTERVALS ON THE AXIS

      CALL PLTINI (X (1), X (2), XSTART, XEND, XINT, IXEXP, IXTIC)
      CALL PLTINI (Y (1), Y (2), YSTART, YEND, YINT, IYEXP, IYTIC)

C  SET ALL THE BUFFER PARAMETERS

      BUFF (1) = 4.
      BUFF (2) = X (1)
      IF (IXEXP .EQ. 0) THEN
         BUFF (3) = XSTART
         BUFF (5) = XINT
      ELSE
         BUFF (3) = XSTART *  (10. ** DBLE(IXEXP))
         BUFF (5) = XINT *  (10. ** DBLE(IXEXP))
      ENDIF
      BUFF (4) = X (2)
      BUFF (6) = 1.
      BUFF (7) = Y (1)
      IF (IYEXP .EQ. 0) THEN
         BUFF (8) = YSTART
         BUFF (10) = YINT
      ELSE
         BUFF (8) = YSTART *  (10. ** DBLE(IYEXP))
         BUFF (10) = YINT *  (10. ** DBLE(IYEXP))
      ENDIF
      BUFF (9) = Y (2)
      BUFF (11) = 1

C  FORCE THE CORRECT AXIS SETUP

      CALL PLTSTG (11, BUFF)

C  PLOT THE AXIS

      CALL PLTGPH (X, Y,  - 2, 'X', ' ', 'Y', ' ')

C  PUT THE CLIPPING RECTANGLE RIGHT AT THE AXIS

      CALL MP2PT (1, X (1), Y (1), X01, Y01, MASK)
      CALL MP2PT (1, X (2), Y (2), X02, Y02, MASK)
      CALL MPVIEW (X01, X02, Y01, Y02)
      CALL MPORT2 (X (1), X (2), Y (1), Y (2))

      RETURN

      END
