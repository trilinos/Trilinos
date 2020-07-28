C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PAINTF (CNTR0, CNTR1, VARNP, NLNKF, LINKF1,
     &   XN, YN, ZN, XF, YF, ZF)
C=======================================================================

C   --*** PAINTF *** (DETOUR) Paint face contour for a face
C   --   Modified by John H. Glick - 10/26/88
C   --   Written by Amy Gilkey - revised 10/22/87
C   --
C   --PAINTF paints a contour section for a face.  The contour algorithm
C   --assumes that the elements do not contain internal nodes.
C   --The element interpolation field is approximated by logically drawing
C   --lines from each node to the element center and, thusly, dividing the
C   --element into triangles.  Contour sections are then drawn by connecting
C   --the intersection points of the sub-element edges and the contour
C   --plane.
C   --
C   --Parameters:
C   --   CNTR0, CNTR1 - the contour values delimiting the contour section
C   --   VARNP - IN - the contour function values
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF1 - IN - the connectivity for the face
C   --   XN, YN, ZN - IN - the nodal coordinates
C   --   XF, YF, ZF - IN - the face center coordinates
C   --
C   --Common Variables:
C   --   Uses IS3DIM of /D3NUMS/

      PARAMETER (KDONE=0, KERR1=-1, KERR2=-2, KERR3=-3,
     &   KCORNR=1, KISIDE=2, KMIDDL=3, KOSIDE=4, KFINDS=5)

      COMMON /D3NUMS/ IS3DIM, NNPSUR, NUMNPF, LLNSET
      LOGICAL IS3DIM

      REAL VARNP(*)
      INTEGER LINKF1(NLNKF)
      REAL XN(*), YN(*), ZN(*)

      LOGICAL INTERP_BL
      LOGICAL INTERC, INTEQC
      LOGICAL ISDONE(8)
      REAL XPTS(24), YPTS(24)

      ICNVS (I) = MOD (I+NNPF-1, NNPF) + 1
C      --ICNVS converts out of range side numbers (1..NLNKF)

      INTERC (CN, C0, C1) = (C0 .LT. CN) .AND. (CN .LT. C1)
C      --INTERC returns true iff corner is within two contours
      INTEQC (CN, CNN, C0, C1) = (CN .EQ. CNN)
     &   .AND. ((CN .EQ. C0) .OR. (CN .EQ. C1))
C      --INTEQC returns true iff corners are equal to each other and a contour

      if (nlnkf .eq. 4) then
C -- Most common case for 3D
C -- Compute center function value and coordinates

        FM = (varnp(linkf1(1)) + varnp(linkf1(2)) + varnp(linkf1(3)) +
     *    varnp(linkf1(4))) / 4.0
         XM = XF
         YM = YF
         NTRI = NLNKF
         NNPF = NLNKF
         MAXPT = 12

       else IF (NLNKF .EQ. 3) THEN

C      --Special case - triangle element

         N = LINKF1(2)
         FM = VARNP(N)
         XM = XN(N)
         YM = YN(N)
         IFL = 0
         NTRI = 1
         NNPF = NLNKF

      ELSE IF ((.NOT. IS3DIM) .AND. (NLNKF .EQ. 8)) THEN

C      --Compute center function value and coordinates

         FM = -0.25 * (varnp(linkf1(1)) + varnp(linkf1(3)) +
     *      varnp(linkf1(5)) + varnp(linkf1(7))) + 0.5 *
     *     (varnp(linkf1(2)) + varnp(linkf1(4)) +
     *      varnp(linkf1(6)) + varnp(linkf1(8)))
         FM = FM
         XM = XF
         YM = YF
         NTRI = NLNKF
         NNPF = NLNKF
         MAXPT = 24

      ELSE IF ((.NOT. IS3DIM) .AND. (NLNKF .EQ. 9)) THEN

C      --Compute center function value and coordinates

         FM = VARNP(LINKF1(9))
         XM = XN(LINKF1(9))
         YM = YN(LINKF1(9))
         NTRI = NLNKF - 1
         NNPF = NLNKF - 1
         MAXPT = 24

      ELSE

C      --Compute center function value and coordinates

         FM = 0.0
         DO 120 K = 1, NLNKF
            FM = FM + VARNP(LINKF1(K))
  120    CONTINUE
         FM = FM / NLNKF
         XM = XF
         YM = YF
         NTRI = NLNKF
         NNPF = NLNKF
         MAXPT = 12
      END IF

C   --Special check to paint the triangle formed by two adjacent nodes
C   --equal to each other, the middle value, and CNTR0

      IF ((NTRI .GT. 1) .AND. (CNTR0 .EQ. FM)) THEN
         IFL = 0
         FLAST = VARNP(LINKF1(NNPF))
         DO 130 K = 1, NNPF
            FTHIS = VARNP(LINKF1(K))
            IF ((FLAST .EQ. FTHIS) .AND. (FTHIS .EQ. FM)) IFL = K
            FLAST = FTHIS
  130    CONTINUE

         IF (IFL .GT. 0) THEN

            XPTS(1) = XM
            YPTS(1) = YM
            N1 = LINKF1(IFL)
            XPTS(2) = XN(N1)
            YPTS(2) = YN(N1)
            N = LINKF1(ICNVS(IFL-1))
            XPTS(3) = XN(N)
            YPTS(3) = YN(N)
            CALL MPD2PG (3, XPTS, YPTS, 'S')
         END IF
      END IF

      CALL INILOG (NTRI, .FALSE., ISDONE)

      NPTS = 0
      DO 150 ISIDE = 1, NTRI

C      --Check each side of the face, looking for a corner in the
C      --contour area or a dividing line in the contour area

         IF (ISDONE(ISIDE)) GOTO 150

         IS = ISIDE
         N = LINKF1(IS)
         N1 = LINKF1(ICNVS(IS+1))
         ISTART = IS

         ISTATE = KDONE
         IF (INTERC (VARNP(N), CNTR0, CNTR1)) THEN
            ICNTR = -1
            ISTATE = KCORNR
         ELSE IF ((INTERP_BL (CNTR0, VARNP(N), VARNP(N1), PSI))
     &      .AND. (VARNP(N1) .GT. CNTR0)) THEN
            ICNTR = 0
            CNTR = CNTR0
            ISTATE = KISIDE
         END IF

  140    CONTINUE
         IF ((ISTATE .GT. KDONE) .AND. (NPTS .LT. MAXPT)) THEN

            IF (ISTATE .EQ. KCORNR) THEN

C            --Corner is within the contour area, goto side-1

               NPTS = NPTS + 1
               XPTS(NPTS) = XN(N)
               YPTS(NPTS) = YN(N)
               ISDONE(IS) = .TRUE.
               IS = ICNVS(IS-1)
               ISDONE(IS) = .TRUE.
               N = LINKF1(IS)
               N1 = LINKF1(ICNVS(IS+1))
               ISTATE = KFINDS
               IF (IS .EQ. ISTART) ISTATE = KDONE

            ELSE IF (ISTATE .EQ. KISIDE) THEN

C            --Side crossing starts a contour area delimiter, find
C            --direction and follow contour line

               NPTS = NPTS + 1
               XPTS(NPTS) = XN(N) * (1.-PSI) + XN(N1) * PSI
               YPTS(NPTS) = YN(N) * (1.-PSI) + YN(N1) * PSI
               ISDONE(IS) = .TRUE.

C            --Find the direction to follow this contour line through
C            --the face center

               ISTATE = KERR3
               IF (INTERP_BL (CNTR, VARNP(N1), FM, PSI)) THEN
                  IDIR = +1
                  I = N1
                  ISTATE = KMIDDL
               ELSE IF (INTERP_BL (CNTR, VARNP(N), FM, PSI)) THEN
                  IDIR = -1
                  I = N
                  ISTATE = KMIDDL
               END IF

            ELSE IF (ISTATE .EQ. KMIDDL) THEN

C            --Crossing through the middle, continue following line
C            --to the side or through the middle

               NPTS = NPTS + 1
               XPTS(NPTS) = XN(I) * (1.-PSI) + XM * PSI
               YPTS(NPTS) = YN(I) * (1.-PSI) + YM * PSI
               ISDONE(IS) = .TRUE.
               IS = ICNVS(IS+IDIR)
               ISDONE(IS) = .TRUE.
               N = LINKF1(IS)
               N1 = LINKF1(ICNVS(IS+1))

               ISTATE = KERR1
               IF (IDIR .EQ. +1) THEN
                  I = N1
               ELSE
                  I = N
               END IF
               IF (INTERP_BL (CNTR, VARNP(N), VARNP(N1), PSI)) THEN
                  ISTATE = KOSIDE
C               --If the corner is exactly on the contour range, check
C               --if the line is moving away from the corner; if so, go
C               --through the middle to prevent cycling on this corner
                  IF (INTERP_BL (CNTR, VARNP(I), FM, PSIX)) THEN
                     IF (((IDIR .EQ. +1) .AND. (VARNP(N) .EQ. CNTR))
     &                  .OR.
     &                  ((IDIR .EQ. -1) .AND. (VARNP(N1) .EQ. CNTR)))
     &                  THEN
                        PSI = PSIX
                        ISTATE = KMIDDL
                     END IF
                  END IF
               ELSE IF (INTERP_BL (CNTR, VARNP(I), FM, PSI)) THEN
                  ISTATE = KMIDDL
               END IF

            ELSE IF (ISTATE .EQ. KOSIDE) THEN

C            --Contour line is intersecting the side, ending the line;
C            --search for enclosing corners/lines

               NPTS = NPTS + 1
               XPTS(NPTS) = XN(N) * (1.-PSI) + XN(N1) * PSI
               YPTS(NPTS) = YN(N) * (1.-PSI) + YN(N1) * PSI
               ISDONE(IS) = .TRUE.
               IF (VARNP(N) .EQ. CNTR) ISDONE(ICNVS(IS-1)) = .TRUE.
               IF (VARNP(N1) .EQ. CNTR) ISDONE(ICNVS(IS+1)) = .TRUE.

               ISTATE = KFINDS
               IF (IS .EQ. ISTART) ISTATE = KDONE
               IF (VARNP(N) .EQ. CNTR) THEN
                  IF ((IDIR .EQ. -1) .OR.
     &               ((ICNTR .EQ. 0) .AND. (IDIR .EQ. +1))) THEN
C                  --If the line exits at the corner, and enclosing area,
C                  --adjust side to be adjacent side
                     IS = ICNVS(IS-1)
                     ISDONE(IS) = .TRUE.
                     N = LINKF1(IS)
                     N1 = LINKF1(ICNVS(IS+1))
                     IF (IS .EQ. ISTART) ISTATE = KDONE
                  END IF
               END IF
            END IF

            IF (ISTATE .EQ. KFINDS) THEN

C            --Search for corner or first side crossing; the side can be
C            --a CNTR0 crossing if the last crossing was a corner or a
C            --CNTR1 crossing; it can be a CNTR1 crossing if the last
C            --crossing was a corner or a CNTR0 crossing

               ISTATE = KERR2
               IF (INTERC (VARNP(N), CNTR0, CNTR1) .OR.
     &            INTEQC (VARNP(N), VARNP(N1), CNTR0, CNTR1)) THEN
                  ICNTR = MIN (-1, ICNTR-1)
                  ISTATE = KCORNR
               ELSE IF ((ICNTR .EQ. -2) .AND.
     &            INTERP_BL (CNTR0, VARNP(N), VARNP(N1), PSI)) THEN
                  ICNTR = 0
                  CNTR = CNTR0
                  ISTATE = KISIDE
               ELSE IF ((ICNTR .NE. 1) .AND.
     &            INTERP_BL (CNTR1, VARNP(N), VARNP(N1), PSI)) THEN
                  ICNTR = 1
                  CNTR = CNTR1
                  ISTATE = KISIDE
               ELSE IF ((ICNTR .NE. 0) .AND.
     &            INTERP_BL (CNTR0, VARNP(N), VARNP(N1), PSI)) THEN
                  ICNTR = 0
                  CNTR = CNTR0
                  ISTATE = KISIDE
               END IF
            END IF

            GOTO 140
         END IF

         IF (ISTATE .NE. KDONE) THEN
C         --State will be <0 if error occurred, >0 if too many points
            CALL PLTFLU
            WRITE (*, 10000) ISTATE, IS, ICNTR, IDIR,
     &         (VARNP(LINKF1(I)), I=1,NLNKF), CNTR0, CNTR1
10000        FORMAT (' PAINT ERROR - please contact',
     &         ' SEACAS@sandia.gov with this information:',
     &         /, 3X, 3I3, SP, I3,
     &         2X, SP, 4 (1X, E9.2), 2X, 2 (1X, E9.2), /)
         END IF

C      --Paint enclosed area

         IF (NPTS .GT. 2) THEN
            CALL MPD2PG (NPTS, XPTS, YPTS, 'S')
            NPTS = 0
         END IF
  150 CONTINUE

      RETURN
      END
