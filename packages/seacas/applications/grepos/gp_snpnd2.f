C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      subroutine snpnd2(numnp, ndim, x, y, vnorm,
     *         neessm, nnessm, ltnesm, toler, delmax)
C=======================================================================

      REAL X(*), Y(*), VNORM(3,*)
      DIMENSION LTNESM(2,*)

      REAL SAVEIT(2)
      LOGICAL FOUND

C ... Quiet the compiler...
      saveit(1) = 0.0
      saveit(2) = 0.0

      glotol = 0.0
      notmat = 0
      matin  = 0
      mattol = 0
      numnod = 0
      dmax = 0.0
      svdel = 0.0

      do 130 inod = 1, numnp
        pmin = 1.0e38
        svmax = -1.0e38
        found = .false.
        if (vnorm(1,inod) .ne. 0.0 .or. vnorm(2,inod) .ne. 0.0) then

          numnod = numnod + 1

          X0 = X(inod)
          Y0 = Y(inod)

          AI = VNORM(1, inod)
          BJ = VNORM(2, inod)
C ... Node movement (delta) = (xnew-X0)**2 + (ynew-Y0)**2
C                           = (x0+t*ai-x0)**2 + (y0+t*bj-y0)**2
C                           = t**2 * (ai**2 + bj**2)
C    Want delta < delmax  ==> t**2 * (ai**2 + bj**2) < delmax**2
C                         ==> t**2 < delmax**2 / (ai**2 + bj**2) = tmax

          tmax = delmax**2 / (ai**2 + bj**2)

          do 110 iseg = 1, neessm
            XI = x(LTNESM(1,ISEG))
            YI = y(LTNESM(1,ISEG))

            XJ = x(LTNESM(2,ISEG))
            YJ = y(LTNESM(2,ISEG))

C ... If denom == 0, then node normal is parallel to plane
            denom = (yj-yi)*ai - (xj-xi)*bj
            if (denom .ne. 0.0) then
              T = ((xj-xi)*(y0-yi) - (yj-yi)*(x0-xi))/denom
              S = (    ai *(y0-yi) -     bj *(x0-xi))/denom

              if (t**2 .le. tmax .and.
     *          0.0 .le. S .and. S .le. 1.0) then
C ... If we made it this far, then the intersection point is inside the
C     face. Move the node to the intersection point and get the next node
                X(INOD) = X0 + T * AI
                Y(INOD) = Y0 + T * BJ

                delta = t**2 * (ai**2 + bj**2)
                dmax = max(delta, dmax)
                matin = matin + 1

                go to 120
              else if (t**2 .le. tmax .and.
     *          0.0-toler .le. S .and. S .le. 1.0+toler) then
C ... If we made it this far, then the intersection point is outside the
C     face, but inside the tolerance range.
C     Save the intersection point in case no intersections within the
C     face are found.
                if (S .gt. 1.0) S = 1.0 - S
C ...           -toler < S < 0.0
                if (toler .gt. svmax) then
                  found = .TRUE.
                  svmax = toler
                  svdel = t**2 * (ai**2 + bj**2)
                  SAVEIT(1) = X0 + T * AI
                  SAVEIT(2) = Y0 + T * BJ
                end if
              else
C ... The node is outside the tolerance,
C     for this face/node combination. Save the minimum for all face/node comb
                if (S .lt. 0.0) then
                  S = -S
                else if (S .gt. 1.0) then
                  S = S - 1.0
                endif

                pmin = min(S, pmin)
              end if
            end if
 110      continue
          if (found) then
            x(inod) = saveit(1)
            y(inod) = saveit(2)
            dmax = max(svdel, dmax)
            mattol = mattol + 1
          else
            write (*,10000) inod, pmin
            glotol = Max(glotol, pmin)
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

