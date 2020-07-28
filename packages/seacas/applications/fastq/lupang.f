C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE LUPANG (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN,
     &   NLOOP, ANGLE, LNODES, NSTART, LLL, XMIN, XMAX, YMIN, YMAX,
     &   ZMIN, ZMAX, DEV1, KREG, ERR)
C***********************************************************************

C  SUROUTINE LUPANG = CALCULATES THE NEW ANGLES FOR ALL NODES IN A LOOP

C***********************************************************************

      DIMENSION XN (MXND), YN (MXND), ZN(MXND)
      DIMENSION LXN(4, MXND), NXL(2, 3*MXND)
      DIMENSION LXK(4, MXND), KXL(2, 3*MXND)
      DIMENSION ANGLE (MXND), LNODES (MLN, MXND)

      LOGICAL ERR

      CHARACTER*3 DEV1

      ERR = .FALSE.

C  LOOP AROUND THE INTERIOR PERIMETER CALCULATING THE NEW
C  ANGLES

      N1 = NSTART
      KOUNT = 0
  100 CONTINUE
      N0 = LNODES (2, N1)
      N2 = LNODES (3, N1)
      CALL GETANG (MXND, MLN, XN, YN, LNODES, LXK, KXL, NXL, LXN,
     &   N0, N1, N2, ANGLE (N1), ERR)
      IF (ERR) THEN
         CALL MESAGE (' ** ERROR IN LUPANG ** ')
         GOTO 110
      ENDIF
      N1 = N2
      IF (N1 .EQ. NSTART) GOTO 110
      KOUNT = KOUNT+1
      IF (KOUNT .GT. NLOOP) THEN
         CALL MESAGE (' ** ERROR IN LUPANG ** ')
         ERR = .TRUE.
         GOTO 110
      ENDIF
      GOTO 100

  110 CONTINUE
      RETURN

      END
