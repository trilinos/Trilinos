C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      subroutine movnod(numnp, ndim, x, y, z,
     *  vnorm, neesss, nnesss, ltness,
     *  plane, neessm, nnessm, ltnesm,
     *  toler, delmax, index, vector, gap)
C=======================================================================

      REAL X(*), Y(*), Z(*), VNORM(3,*), PLANE(4,*)
      DIMENSION LTNESS(4,*), LTNESM(4,*)
      INTEGER INDEX(*)
      REAL VECTOR(3)
      LOGICAL FOUND

      REAL PI(3)

      dmax = 0.0
      dmin = 1.0e15
      delmx2 = delmax**2
      match = 0

      do 130 inod = 1, numnp
        found = .FALSE.
        smin = 1.0e15
        smax = 0.0
        armax = -1.0e38
        if (vnorm(1,inod) .ne. 0.0 .or. vnorm(2,inod) .ne. 0.0 .or.
     *    vnorm(3,inod) .ne. 0.0) then

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

              if (t .lt. 0.0) go to 100
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
C     face. Save the minimum distance found so far.
                  dmin = min(delta2, dmin)
                  dmax = max(delta2, dmax)
                  match = match + 1

                  go to 120
                else if (area1 .ge. -toler .and. area2 .ge. -toler .and.
     *            area3 .ge. -toler .and. area4 .ge. -toler) then
C ... If we made it this far, then the intersection point is outside the
C     face, but inside the toler. Save the minimum distance found so far.
                  found = .true.
                  smin = min(delta2, smin)
                  smax = max(delta2, smax)
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
C ... The node did not intersect any face, but it did find one within
C     the tolerance. Use the closest one
            dmin = smin
          end if
        end if
 120    continue
 130  continue

C ... Update the node positions based on the minimum distance found
C     and the specified vector.
      if (match .gt. 0) then
        AI = vector(1)
        BJ = vector(2)
        CK = vector(3)

        dmin = sqrt(dmin)
        dmin = dmin - gap

        do 140 inod=1, numnp
          if (index(inod) .eq. 1) then
            X0 = X(inod)
            Y0 = Y(inod)
            Z0 = Z(inod)

C ... Update the nodes position (Currently, assumes in vector direction)
            X(inod) = X0 + dmin * AI
            Y(inod) = Y0 + dmin * BJ
            Z(inod) = Z0 + dmin * CK
          end if
 140    continue
        write (*, 10020) dmin
      else
        write (*,*) 'No node movement.'
      end if

10020 format(/,'Node movement = ',1pe11.4)
      return
      end

