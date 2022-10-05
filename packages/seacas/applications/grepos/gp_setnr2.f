C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SETNR2 (ISPLANE, NUMNP, X, Y, IDIM, VNORM, NRMTYP,
     *  PRJVEC, NEESS, NNESS, LTNESS)
C=======================================================================

      include 'gp_snap.blk'

      LOGICAL ISPLANE
      INTEGER NRMTYP

      REAL X(*), Y(*), VNORM(IDIM,*), PRJVEC(3)
      DIMENSION LTNESS(2,*)

      IF (NRMTYP .EQ. PNORM .OR. ISPLANE) THEN
        DO 20 ISEG = 1, NEESS

            XI = x(LTNESS(1,ISEG))
            YI = y(LTNESS(1,ISEG))

            XJ = x(LTNESS(2,ISEG))
            YJ = y(LTNESS(2,ISEG))

            DX = XI - XJ
            DY = YI - YJ
            RMAG = SQRT ( DX**2 + DY**2)

            AI = -dy / rmag
            BJ =  dx / rmag

          if (isplane) then
            vnorm(1, iseg) = ai
            vnorm(2, iseg) = bj
            vnorm(3, iseg) = 0.0
            di = ai * xi + bj * yi
            dj = ai * xj + bj * yj
            vnorm(4,iseg) = (di + dj) / 2.0
          else
            DO 10 INOD = 1, 2
              NODE = LTNESS(INOD, ISEG)
              VNORM(1, NODE) = VNORM(1, NODE) + AI
              VNORM(2, NODE) = VNORM(2, NODE) + BJ
 10         CONTINUE
          end if
 20     CONTINUE

        if (isplane) return

C ... NORMALIZE PRJVECS
        do 30 i=1, neess
          do 25 j=1,2
            node = ltness(j, i)
            A = VNORM(1,NODE)
            B = VNORM(2,NODE)
            R = A**2 + B**2
            IF (R .GT. 0.0) THEN
              R = SQRT(R)
              VNORM(1,NODE) = A/R
              VNORM(2,NODE) = B/R
            END IF
 25       continue
 30     CONTINUE

      ELSE IF (NRMTYP .EQ. PRAD) THEN
C ... Radial Projection
        do 70 i=1, neess
          do 60 j=1,2
            node = ltness(j, i)
            a = x(node) - prjvec(1)
            b = y(node) - prjvec(2)
            R = sqrt(a**2 + b**2)
            if (r .eq. 0) r = 1.0
            vnorm(1,node) = a/r
            vnorm(2,node) = b/r
 60       continue
 70     continue

      ELSE IF (NRMTYP .EQ. PVECT) THEN
C ... Use user-supplied normal
        R = SQRT(PRJVEC(1)**2 + PRJVEC(2)**2)
        R = SQRT(R)
        PRJVEC(1) = PRJVEC(1)/R
        PRJVEC(2) = PRJVEC(2)/R
        do 50 i=1, neess
          do 40 j=1,2
            node = ltness(j, i)
            vnorm(1,node) = prjvec(1)
            vnorm(2,node) = prjvec(2)
 40       continue
 50     continue

      ELSE IF (NRMTYP .EQ. PEDGE) THEN

      END IF
      RETURN
      END

