C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      subroutine snpnod(numnp, ndim, x, y, z,
     *  vnorm, neesss, nnesss, ltness,
     *  plane, neessm, nnessm, ltnesm,
     *  toler, delmax)
C=======================================================================

      REAL X(*), Y(*), Z(*), VNORM(3,*), PLANE(4,*)
      DIMENSION LTNESS(4,*), LTNESM(4,*)

      REAL PI(3)
      REAL SAVEIT(3)
      LOGICAL FOUND

      glotol = 0.0
      notmat = 0
      matin  = 0
      mattol = 0
      numnod = 0
      dmax = 0.0
      delmx2 = delmax**2
      svdel = 0.0
      saveit(1) = 0.0
      saveit(2) = 0.0
      saveit(3) = 0.0

      do 130 inod = 1, numnp
        armax = -1.0e38
        svmax = -1.0e38
        found = .false.
        if (vnorm(1,inod) .ne. 0.0 .or. vnorm(2,inod) .ne. 0.0 .or.
     *    vnorm(3,inod) .ne. 0.0) then

          numnod = numnod + 1

          X0 = X(inod)
          Y0 = Y(inod)
          Z0 = Z(inod)

          AI = VNORM(1, inod)
          BJ = VNORM(2, inod)
          CK = VNORM(3, inod)

          do 110 ifac = 1, neessm
            A = plane(1,ifac)
            B = plane(2,ifac)
            C = plane(3,ifac)
            D = plane(4,ifac)

C ... If denom == 0, then node normal is parallel to plane
            DENOM = A*AI + B*BJ + C*CK
            if (denom .ne. 0.0) then
              T = -(A*X0 + B*Y0 + C*Z0 - D) / DENOM

C ... Intersection point
              PI(1) = X0 + T * AI
              PI(2) = Y0 + T * BJ
              PI(3) = Z0 + T * CK

              dx = abs(x(inod)-pi(1))
              if (dx .gt. delmax) go to 100
              dy = abs(y(inod)-pi(2))
              if (dy .gt. delmax) go to 100
              dz = abs(z(inod)-pi(3))
              if (dz .gt. delmax) go to 100

              delta2 = dx**2 + dy**2 + dz**2
              if (delta2 .le. delmx2) then
C ... See if intersection point is inside face.
                XI = X(LTNESM(1,IFAC))
                YI = Y(LTNESM(1,IFAC))
                ZI = Z(LTNESM(1,IFAC))

                XJ = X(LTNESM(2,IFAC))
                YJ = Y(LTNESM(2,IFAC))
                ZJ = Z(LTNESM(2,IFAC))

                XK = X(LTNESM(3,IFAC))
                YK = Y(LTNESM(3,IFAC))
                ZK = Z(LTNESM(3,IFAC))

                XL = X(LTNESM(4,IFAC))
                YL = Y(LTNESM(4,IFAC))
                ZL = Z(LTNESM(4,IFAC))

                area1 = trarea(xi, yi, zi, xj, yj, zj,
     *            pi(1), pi(2), pi(3),
     *            plane(1,ifac))

                area2 = trarea(xj, yj, zj, xk, yk, zk,
     *            pi(1), pi(2), pi(3),
     *            plane(1,ifac))

                area3 = trarea(xk, yk, zk, xl, yl, zl,
     *            pi(1), pi(2), pi(3),
     *            plane(1,ifac))

                area4 = trarea(xl, yl, zl, xi, yi, zi,
     *            pi(1), pi(2), pi(3),
     *            plane(1,ifac))

                if (area1 .ge. 0.0 .and. area2 .ge. 0.0 .and.
     *              area3 .ge. 0.0 .and. area4 .ge. 0.0) then
C ... If we made it this far, then the intersection point is inside the
C     face. Move the node to the intersection point and get the next node
                  X(INOD) = PI(1)
                  Y(INOD) = PI(2)
                  Z(INOD) = PI(3)

                  dmax = max(delta2, dmax)
                  matin = matin + 1

                  go to 120
                else if (area1 .ge. -toler .and. area2 .ge. -toler .and.
     *            area3 .ge. -toler .and. area4 .ge. -toler) then
C ... If we made it this far, then the intersection point is slightly
C     outside the face, but is within tolerance. Save this intersection
C     point in case no intersections within the face are found.

C ...             Want largest negative
                  armin = MIN(area1, area2, area3, area4)

C ...             NOTE: armin < 0, want closest to zero.
                  if (armin .gt. svmax) then
                    found = .TRUE.
                    svmax = armin
                    svdel = delta2
                    SAVEIT(1) = PI(1)
                    SAVEIT(2) = PI(2)
                    SAVEIT(3) = PI(3)
                  end if
                else
C ... The node is outside the tolerance, find the most negative of the area*
C     for this face/node combination. Save the maximum (closest to zero)
C     for all face/node combinations.
                  armin = MIN(area1, area2, area3, area4)
                  armax = MAX(armin, armax)
                end if
              end if
            end if
 100        continue
 110      continue
          if (found) then
            x(inod) = saveit(1)
            y(inod) = saveit(2)
            z(inod) = saveit(3)
            dmax = max(svdel, dmax)
            mattol = mattol + 1
          else
            write (*,10000) inod, -armax
            glotol = MAX(glotol, -armax)
            notmat = notmat + 1
          end if
        end if
 120    continue
 130  continue

      if (notmat .gt. 0) then
        write (*,10010) notmat, glotol
      end if
      write (*, 10020) sqrt(dmax), numnod, matin, mattol
10000 format('WARNING - No matching face found for node ',i9,/,
     *       '          Tolerance must be at least ',1pe11.4,
     *       ' to detect a match.')
10010 format(/,'WARNING - ',I9,' nodes were not matched.',/,
     *  'Set tolerance greater than ',1pe11.4,' to match all nodes.')
10020 format(/,'Maximum node movement        = ',1pe11.4,
     *       /,'Number of unique slave nodes = ',I9,
     *       /,'Number of exact matches      = ',I9,
     *       /,'Number of toleranced matches = ',I9 )
      return
      end

      real function trarea(A1, A2, A3,   B1, B2, B3,  C1, C2, C3,
     *  RNORM)
      real sum(3), rnorm(3)

      sum(1) = A2*B3 - A3*B2
      sum(2) = A3*B1 - A1*B3
      sum(3) = A1*B2 - A2*B1

      sum(1) = sum(1) + B2*C3 - B3*C2
      sum(2) = sum(2) + B3*C1 - B1*C3
      sum(3) = sum(3) + B1*C2 - B2*C1

      sum(1) = sum(1) + C2*A3 - C3*A2
      sum(2) = sum(2) + C3*A1 - C1*A3
      sum(3) = sum(3) + C1*A2 - C2*A1

      trarea = sum(1)*rnorm(1) + sum(2)*rnorm(2) + sum(3)*rnorm(3)
      return
      end

