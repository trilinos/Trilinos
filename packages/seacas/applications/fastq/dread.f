C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE DREAD (X, Y, BUTTON, X1, Y1, XX1, YY1, SCALE, CT, ST)
C***********************************************************************

C  SUBROUTINE DREAD = SETS ALL PARAMETERS UP FOR READING FROM A DIGI-PAD

C***********************************************************************

C  SUBROUTINES CALLED:
C     DPREAD   = READS INPUT FROM A DIGI-PAD DIGITIZER

C***********************************************************************

C  VARIABLES USED:
C     X      = THE X LOCATION IN USER COORDINATES
C     Y      = THE Y LOCATION IN USER COORDINATES
C     BUTTON = THE MOUSE BUTTON PUSHED
C     X1     = THE X ZERO SHIFT IN USER COORDINATES
C     Y1     = THE Y ZERO SHIFT IN USER COORDINATES
C     XX1    = THE X ZERO SHIFT IN DIGITIZED COORDINATES
C     YY1    = THE Y ZERO SHIFT IN DIGITIZED COORDINATES
C     SCALE  = THE SCALE FACTOR FROM DIGITIZED TO USER COORDINATES
C     CT     = THE COSINE OF THE ANGLE OF THE DRAWING ON THE PAD
C     ST     = THE SINE OF THE ANGLE OF THE DRAWING ON THE PAD
C     XNEW   = THE NEW DIGITIZED X VALUE BEFORE TRANSFORMATIONS
C     YNEW   = THE NEW DIGITIZED Y VALUE BEFORE TRANSFORMATIONS
C     NCB    = THE NUMBER OF BUTTONS ON THE MOUSE  (BIT-PAD-ONE)
C     DEL    = THE DELTA DISTANCE BETWEEN ACCEPTABLE POINTS  (TALOS)

C***********************************************************************

      CHARACTER*1 BUTTON

      CALL DPREAD (XNEW, YNEW, BUTTON)
      XNEW = XNEW - XX1
      YNEW = YNEW - YY1
      X = X1 + SCALE * ( (CT * XNEW) + (ST * YNEW))
      Y = Y1 + SCALE * ( (-ST * XNEW) + (CT * YNEW))
      RETURN
      END
