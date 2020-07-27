C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE FLMNMX (MXND, MLN, MAXPRM, LINKPR, KPERIM, LNODES,
     &   XN, YN, NLOOP, NODE, XMIN, XMAX, YMIN, YMAX, ERR)
C***********************************************************************

C  SUBROUTINE FLMNMX = SET MIN AND MAX FOR CURRENT FILL BOUNDARY

C***********************************************************************

      DIMENSION LNODES (MLN, MXND), XN (MXND), YN (MXND)
      DIMENSION LINKPR (3, MAXPRM)

      LOGICAL ERR

      KOUNT = 0
      INOW = NODE
      XMIN = XN (NODE)
      XMAX = XN (NODE)
      YMIN = YN (NODE)
      YMAX = YN (NODE)

  100 CONTINUE

      INOW = LNODES (3, INOW)
      IF (INOW .NE. NODE) THEN

         XMIN = MIN (XMIN, XN (INOW))
         YMIN = MIN (YMIN, YN (INOW))
         XMAX = MAX (XMAX, XN (INOW))
         YMAX = MAX (YMAX, YN (INOW))

         KOUNT = KOUNT + 1

         IF (KOUNT .GT. NLOOP) THEN
            CALL MESAGE('PROBLEMS IN FLMNMX WITH LOOP NOT CLOSING')
            ERR = .TRUE.
            GOTO 130
         ENDIF
         GOTO 100
      ENDIF

C  LOOP THROUGH ALL THE REMAINING PERIMETERS CHECKING FOR CROSSINGS

      IPERIM = KPERIM
  110 CONTINUE
      IPERIM = LINKPR (2, IPERIM)
      IF ((IPERIM .EQ. 0) .OR. (IPERIM .EQ. KPERIM)) GOTO 130

      KMAX = LINKPR (3, IPERIM)
      INOW = LINKPR (1, IPERIM)
      KOUNT = 0

  120 CONTINUE
      XMIN = MIN (XMIN, XN (INOW))
      YMIN = MIN (YMIN, YN (INOW))
      XMAX = MAX (XMAX, XN (INOW))
      YMAX = MAX (YMAX, YN (INOW))
      KOUNT = KOUNT + 1
      INOW = LNODES (3, INOW)
      IF (INOW .EQ. LINKPR (1, IPERIM)) GOTO 110

      IF (KOUNT. GT. KMAX + 1) THEN
         CALL MESAGE('PROBLEMS IN FLMNMX WITH LOOP NOT CLOSING')
         ERR = .TRUE.
         GOTO 130
      ENDIF
      GOTO 120

  130 CONTINUE
      RETURN

      END
