C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SETNOR (ISPLANE, NUMNP, X, Y, Z, IDIM, VNORM, NRMTYP,
     *  PRJVEC, NEESS, NNESS, LTNESS)
C=======================================================================

      include 'gp_snap.blk'

      LOGICAL ISPLANE
      INTEGER NRMTYP

      REAL X(*), Y(*), Z(*), VNORM(IDIM,*), PRJVEC(3)
      DIMENSION LTNESS(4,*)

      IF (NRMTYP .EQ. PNORM .OR. ISPLANE) THEN
        DO 20 ISEG = 1, NEESS
          XI = X(LTNESS(1,ISEG))
          YI = Y(LTNESS(1,ISEG))
          ZI = Z(LTNESS(1,ISEG))

          XJ = X(LTNESS(2,ISEG))
          YJ = Y(LTNESS(2,ISEG))
          ZJ = Z(LTNESS(2,ISEG))

          XK = X(LTNESS(3,ISEG))
          YK = Y(LTNESS(3,ISEG))
          ZK = Z(LTNESS(3,ISEG))

          XL = X(LTNESS(4,ISEG))
          YL = Y(LTNESS(4,ISEG))
          ZL = Z(LTNESS(4,ISEG))

          AI =  (YK - YI) * (ZL - ZJ) - (ZK - ZI) * (YL - YJ)
          BJ =  (ZK - ZI) * (XL - XJ) - (XK - XI) * (ZL - ZJ)
          CK =  (XK - XI) * (YL - YJ) - (YK - YI) * (XL - XJ)
          RMAG = SQRT ( AI**2 + BJ**2 + CK**2)

          AI = AI / RMAG
          BJ = BJ / RMAG
          CK = CK / RMAG

          if (isplane) then
            vnorm(1, iseg) = ai
            vnorm(2, iseg) = bj
            vnorm(3, iseg) = ck
            di = ai * xi + bj * yi + ck * zi
            dj = ai * xj + bj * yj + ck * zj
            dk = ai * xk + bj * yk + ck * zk
            dl = ai * xl + bj * yl + ck * zl
            vnorm(4,iseg) = (di + dj + dk + dl) / 4.0
          else
            DO 10 INOD = 1, 4
              NODE = LTNESS(INOD, ISEG)
              VNORM(1, NODE) = VNORM(1, NODE) + AI
              VNORM(2, NODE) = VNORM(2, NODE) + BJ
              VNORM(3, NODE) = VNORM(3, NODE) + CK
 10         CONTINUE
          end if
 20     CONTINUE

        if (isplane) return

C ... NORMALIZE PRJVECS
        do 30 i=1, neess
          do 25 j=1,4
            node = ltness(j, i)
            A = VNORM(1,NODE)
            B = VNORM(2,NODE)
            C = VNORM(3,NODE)
            R = A**2 + B**2 + C**2
            IF (R .GT. 0.0) THEN
              R = SQRT(R)
              VNORM(1,NODE) = A/R
              VNORM(2,NODE) = B/R
              VNORM(3,NODE) = C/R
            END IF
 25       continue
 30     CONTINUE

      ELSE IF (NRMTYP .EQ. PRAD) THEN
C ... Radial Projection
        do 70 i=1, neess
          do 60 j=1,4
            node = ltness(j, i)
            a = x(node) - prjvec(1)
            b = y(node) - prjvec(2)
            c = z(node) - prjvec(3)
            R = sqrt(a**2 + b**2 + c**2)
            if (r .eq. 0) r = 1.0
            vnorm(1,node) = a/r
            vnorm(2,node) = b/r
            vnorm(3,node) = c/r
 60       continue
 70     continue

      ELSE IF (NRMTYP .EQ. PVECT) THEN
C ... Use user-supplied normal
        R = SQRT(PRJVEC(1)**2 + PRJVEC(2)**2 + PRJVEC(3)**2)
        PRJVEC(1) = PRJVEC(1)/R
        PRJVEC(2) = PRJVEC(2)/R
        PRJVEC(3) = PRJVEC(3)/R
        do 50 i=1, neess
          do 40 j=1,4
            node = ltness(j, i)
            vnorm(1,node) = prjvec(1)
            vnorm(2,node) = prjvec(2)
            vnorm(3,node) = prjvec(3)
 40       continue
 50     continue

      ELSE IF (NRMTYP .EQ. PEDGE) THEN

      END IF
      RETURN
      END
