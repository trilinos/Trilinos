C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C     -*- Mode: fortran -*-
C=======================================================================

      SUBROUTINE MATXYZ (NDIM,
     &     MATNS1, MATNS2, NNNPS, IXNNPS, LTNNPS,
     &     NUMNP1, XN1, YN1, ZN1, NUMNP2, XN2, YN2, ZN2,
     &     IX1, IX2, IXNP2, NMATCH, TOLER, CLOSE, IEXPCT)
C=======================================================================

C     --*** MATXYZ *** (GJOIN) Find matching nodes
C     --   Written by Amy Gilkey - revised 03/04/88
C     --
C     --MATXYZ matches nodes from the first and second database by comparing
C     --the coordinates.  Since some slack is needed in the equality check,
C     --the comparison may not work for all meshes.  The nodes to be matched
C     --may be limited to nodes in two nodal point sets.
C     --
C     --Parameters:
C     --   NDIM - IN - the number of coordinates per node
C     --   MATNS1 - IN - the number of the first nodal point set
C     --   MATNS2 - IN - the number of the second nodal point set
C     --   NNNPS - IN - the number of nodes for each set
C     --   IXNNPS - IN - the indices of the nodes for each set
C     --   LTNNPS - IN - the nodes for each set
C     --   LNPS2 - IN - the length of the nodal point set for the second database
C     --   NPS2 - IN - the nodes in the nodal point set for the second database
C     --   NUMNP1 - IN - the number of nodes in the first database
C     --   XN1, YN1, ZN1 - IN - the coordinates in the first database
C     --   NUMNP2 - IN - the number of nodes in the second database
C     --   XN2, YN2, ZN2 - IN - the coordinates in the second database
C     --   IX1 - SCRATCH - size = NUMNP1
C     --   IX2 - SCRATCH - size = NUMNP2
C     --   IXNP2 - IN/OUT - the index of the matching first database node;
C     --      negative if no match; reset if NMATCH = 0
C     --   NMATCH - IN/OUT - the number of matching nodes
C     --   CLOSE - IN/OUT - true if match closest node, false if match any
C     --      node within tolerance.
C     --   IEXPCT - IN - number of expected nodal point matches.  Abort
C     --      if actual matches not equal.

      include 'gj_filnum.blk'
      INTEGER NNNPS(*), IXNNPS(*), LTNNPS(*)
      REAL XN1(*), YN1(*), ZN1(*)
      REAL XN2(*), YN2(*), ZN2(*)
      INTEGER IX1(*), IX2(*)
      INTEGER IXNP2(*)
      LOGICAL CLOSE

      PARAMETER (MXNAM = 1)
      DIMENSION KV(MXNAM), RV(MXNAM), IVAL(MXNAM)
      CHARACTER*8  CV(MXNAM)

      CHARACTER*80 STRING
      LOGICAL INIT
      LOGICAL BYSET
      LOGICAL OK

      BYSET = (MATNS1 .GT. 0) .AND. (MATNS2 .GT. 0)
      INIT = (NMATCH .LE. 0)
      NSVMAT = NMATCH

C     --Index nodes to match in the nodal point sets

      IF (BYSET) THEN
         IN1 = 0
         IX0 = IXNNPS(MATNS1) - 1
         DO 100 N = 1, NNNPS(MATNS1)
            IN1 = IN1 + 1
            IX1(IN1) = LTNNPS(IX0+N)
 100     CONTINUE

         IN2 = 0
         IX0 = IXNNPS(MATNS2) - 1
         DO 110 N = 1, NNNPS(MATNS2)
            IN2 = IN2 + 1
            IX2(IN2) = LTNNPS(IX0+N)
 110     CONTINUE
      END IF

      time0 = 0.0
      time1 = 0.0
      time2 = 0.0
C     --Find the limits of the overlapping area of the two databases

      IF (BYSET) THEN
         CALL MINMXS (IN1, IX1, XN1, X1MIN, X1MAX)
         CALL MINMXS (IN2, IX2, XN2, X2MIN, X2MAX)
         CALL MINMXS (IN1, IX1, YN1, Y1MIN, Y1MAX)
         CALL MINMXS (IN2, IX2, YN2, Y2MIN, Y2MAX)
         IF (NDIM .GE. 3) THEN
            CALL MINMXS (IN1, IX1, ZN1, Z1MIN, Z1MAX)
            CALL MINMXS (IN2, IX2, ZN2, Z2MIN, Z2MAX)
         END IF
      ELSE
         CALL MINMAX (NUMNP1, XN1, X1MIN, X1MAX)
         CALL MINMAX (NUMNP2, XN2, X2MIN, X2MAX)
         CALL MINMAX (NUMNP1, YN1, Y1MIN, Y1MAX)
         CALL MINMAX (NUMNP2, YN2, Y2MIN, Y2MAX)
         IF (NDIM .GE. 3) THEN
            CALL MINMAX (NUMNP1, ZN1, Z1MIN, Z1MAX)
            CALL MINMAX (NUMNP2, ZN2, Z2MIN, Z2MAX)
         END IF
      END IF
      XMIN = MAX (X1MIN, X2MIN)
      XMAX = MIN (X1MAX, X2MAX)
      YMIN = MAX (Y1MIN, Y2MIN)
      YMAX = MIN (Y1MAX, Y2MAX)
      IF (NDIM .GE. 3) THEN
         ZMIN = MAX (Z1MIN, Z2MIN)
         ZMAX = MIN (Z1MAX, Z2MAX)
      ELSE
         ZMIN = 0.0
         ZMAX = 0.0
      END IF

      DELTAX = XMAX - XMIN
      DELTAY = YMAX - YMIN
      DELTAZ = ZMAX - ZMIN

      EPS = (DELTAX + DELTAY + DELTAZ) / 1E3
      IF (EPS .LT. 0.0) THEN
         CALL PRTERR ('ERROR', 'Nodes do not overlap')
         GOTO 180
      END IF

      IF (TOLER .GE. 0.0) THEN
         EPS = TOLER
      ELSE
         WRITE (*, 10040) EPS
10040    FORMAT (/' Default tolerance = ',1PE10.3)
         CALL FREFLD (0, 0,
     *        'Enter new value for tolerance (<ret> for default): ',
     *        MXNAM, IOS, NF, KV, CV, IVAL, RV)
         IF (IOS .NE. 0) RV(1) = 0.0
         IF (RV(1) .NE. 0.0) EPS = RV(1)
         CALL OUTLOG (KLOG, 1, 1, CV, IVAL, EPS)
      ENDIF
      OK = .TRUE.
      OK = OK .AND.
     &     ((MIN (X1MAX, X2MAX) - MAX (X1MIN, X2MIN)) .GT. -EPS)
      OK = OK .AND.
     &     ((MIN (Y1MAX, Y2MAX) - MAX (Y1MIN, Y2MIN)) .GT. -EPS)
      IF (NDIM .EQ. 3) THEN
         OK = OK .AND.
     &        ((MIN (Z1MAX, Z2MAX) - MAX (Z1MIN, Z2MIN)) .GT. -EPS)
      END IF
      IF (.NOT. OK) THEN
         CALL PRTERR ('ERROR', 'Nodes do not overlap')
         GOTO 180
      END IF

      XMIN = XMIN - EPS
      XMAX = XMAX + EPS
      YMIN = YMIN - EPS
      YMAX = YMAX + EPS
      ZMIN = ZMIN - EPS
      ZMAX = ZMAX + EPS

C     --Index the nodes within the overlap area

      Z3D = 0.0
      IF (.NOT. BYSET) THEN
         IN1 = 0
         Z3D = 0.0
         DO 120 INP = 1, NUMNP1
            IF ((XN1(INP) .GE. XMIN) .AND. (XN1(INP) .LE. XMAX)) THEN
               IF ((YN1(INP) .GE. YMIN) .AND. (YN1(INP) .LE. YMAX)) THEN
                  IF (NDIM .GE. 3) Z3D = ZN1(INP)
                  IF ((Z3D .GE. ZMIN) .AND. (Z3D .LE. ZMAX)) THEN
                     IN1 = IN1 + 1
                     IX1(IN1) = INP

                  END IF
               END IF
            END IF
 120     CONTINUE

         IN2 = 0
         IF (NDIM .LT. 3) Z3D = 0.0
         DO 130 INP = 1, NUMNP2
            IF ((XN2(INP) .GE. XMIN) .AND. (XN2(INP) .LE. XMAX)) THEN
               IF ((YN2(INP) .GE. YMIN) .AND. (YN2(INP) .LE. YMAX)) THEN
                  IF (NDIM .GE. 3) Z3D = ZN2(INP)
                  IF ((Z3D .GE. ZMIN) .AND. (Z3D .LE. ZMAX)) THEN
                     IN2 = IN2 + 1
                     IX2(IN2) = INP
                  END IF
               END IF
            END IF
 130     CONTINUE
      ELSE

C     ... Equivalencing being done by nodesets.
C     ... Index nodes in nodeset that are within overlap area.
         IN1 = 0
         IX0 = IXNNPS(MATNS1) - 1
         DO 135 N = 1, NNNPS(MATNS1)
            INP = LTNNPS(IX0+N)
            IF ((XN1(INP) .GE. XMIN) .AND. (XN1(INP) .LE. XMAX)) THEN
               IF ((YN1(INP) .GE. YMIN) .AND. (YN1(INP) .LE. YMAX)) THEN
                  IF (NDIM .GE. 3) Z3D = ZN1(INP)
                  IF ((Z3D .GE. ZMIN) .AND. (Z3D .LE. ZMAX)) THEN
                     IN1 = IN1 + 1
                     IX1(IN1) = INP
                  END IF
               END IF
            END IF
 135     CONTINUE

         IN2 = 0
         IX0 = IXNNPS(MATNS2) - 1
         DO 136 N = 1, NNNPS(MATNS2)
            INP = LTNNPS(IX0+N)
            IF ((XN2(INP) .GE. XMIN) .AND. (XN2(INP) .LE. XMAX)) THEN
               IF ((YN2(INP) .GE. YMIN) .AND. (YN2(INP) .LE. YMAX)) THEN
                  IF (NDIM .GE. 3) Z3D = ZN2(INP)
                  IF ((Z3D .GE. ZMIN) .AND. (Z3D .LE. ZMAX)) THEN
                     IN2 = IN2 + 1
                     IX2(IN2) = INP
                  END IF
               END IF
            END IF
 136     CONTINUE
      END IF

C     --Blank out the IXNP2 array

      IF (INIT) THEN
         DO 140 INP = 1, NUMNP2
            IXNP2(INP) = -999
 140     CONTINUE
      ELSE
         DO 150 INP = 1, NUMNP2
            IF (IXNP2(INP) .GT. 0) THEN
               I = LOCINT (INP, IN2, IX2)
               IF (I .GT. 0) THEN
                  IX2(I) = IX2(IN2)
                  IN2 = IN2 - 1
               END IF
               I = LOCINT (IXNP2(INP), IN1, IX1)
               IF (I .GT. 0) THEN
                  IX1(I) = IX1(IN1)
                  IN1 = IN1 - 1
               END IF
            END IF
 150     CONTINUE
      END IF

C     -- Find all matching nodes by comparing coordinates of nodes in overlap
C     area
      DELMAX = MAX(DELTAX, DELTAY, DELTAZ)

      call excpus(time0)
      imat = 0
      if (DELMAX .EQ. DELTAX) THEN
        imat = 1
        call indexn(xn1, 1, 1, ix1, in1, .FALSE.)
        call indexn(xn2, 1, 1, ix2, in2, .FALSE.)
        call prterr('CMDSPEC', 'Entering Sorting Phase, Sort on X')
      else if (delmax .eq. deltay) then
        imat = 2
        call indexn(yn1, 1, 1, ix1, in1, .FALSE.)
        call indexn(yn2, 1, 1, ix2, in2, .FALSE.)
        call prterr('CMDSPEC', 'Entering Sorting Phase, Sort on Y')
      else if (ndim .ge. 3 .and. delmax .eq. deltaz) then
        imat = 3
        call indexn(zn1, 1, 1, ix1, in1, .FALSE.)
        call indexn(zn2, 1, 1, ix2, in2, .FALSE.)
        call prterr('CMDSPEC', 'Entering Sorting Phase, Sort on Z')
      end if
      call prterr('CMDSPEC', 'Entering Comparison Phase')

      CALL EXCPUS(time1)
      DISMIN =  1.0E38
      DISMAX = -1.0E38
      NCOMP  = 0
      nout   = 0
      WRITE (*,'(A)') ' '
      WRITE (*,'(A)') ' '

      Z3D = 0.0
      Z = 0.0
      IN1SV = IN1
      IN2SV = IN2

      I2BEG = 1
      I1BEG = 1

      if (in1sv .lt. in2sv) then

         DO 170 I1 = 1, IN1SV
            INP1 = IX1(I1)
            X = XN1(INP1)
            Y = YN1(INP1)
            IF (NDIM .GE. 3) Z = ZN1(INP1)
            DMIN = 1.0E38
            NDMIN = 0
            I2BEGS = I2BEG
            DO 160 I2 = I2BEGS, IN2SV
               NCOMP = NCOMP + 1
               nout  = nout  + 1
               IF ( nout .ge. 100000) THEN
                  NCMPX = (IN1SV - I1) * IN2 + (IN2 - I2)
                  WRITE (*, 10020) NCOMP, NCMPX
                  nout = 0
               END IF
               INP2 = IX2(I2)
               if (inp2 .le. 0) goto 160
               IF (NDIM .GE. 3) Z3D = ZN2(INP2)

               if ((imat .eq. 1 .and. (x-eps .gt. xn2(inp2))) .or.
     *           (imat .eq. 2 .and. (y-eps .gt. yn2(inp2))) .or.
     *           (imat .eq. 3 .and. (z-eps .gt. z3d))) then
                 i2beg = i2
                 goto 160
               end if

C     ... Since we are sorted on coordinat X|Y|Z,
C     if set 2 X|Y|Z greater than set 1 X|Y|Z+eps, go to next X1|Y1|Z1 coord.
               if ((imat .eq. 1 .and. xn2(inp2)-eps .gt. x) .or.
     *           (imat .eq. 2 .and. yn2(inp2)-eps .gt. y) .or.
     *           (imat .eq. 3 .and. z3d      -eps .gt. z)) goto 165

               DIS = MAX (ABS (XN2(INP2) - X), ABS(YN2(INP2) - Y),
     *              ABS (Z3D - Z) )
               IF ( (DIS .LE. EPS .AND. .NOT. CLOSE) .OR.
     $              DIS .EQ. 0.0) THEN
                  DISMAX = MAX(DISMAX, DIS)
                  NMATCH = NMATCH + 1
                  IXNP2(INP2) = INP1
                  IX2(I2) = -IX2(I2)
                  IN2 = IN2 - 1
                  IX1(I1) = 0
                  IN1 = IN1 - 1
                  GOTO 170
               ELSE IF (DIS .LE. EPS .AND. CLOSE) THEN
                  IF (DIS .LT. DMIN) THEN
                     DMIN = DIS
                     NDMIN = I2
                  END IF
               ELSE
                  DISMIN = MIN(DISMIN, DIS)
               END IF
 160        CONTINUE
 165        CONTINUE
            IF (CLOSE) THEN
               IF (DMIN .LE. EPS) THEN
                  DISMAX = MAX(DISMAX, DMIN)
                  NMATCH = NMATCH + 1
                  INP2 = IX2(NDMIN)
                  IXNP2(INP2) = INP1
                  IX2(NDMIN) = -IX2(NDMIN)
                  IN2 = IN2 - 1
                  IX1(I1) = 0
                  IN1 = IN1 - 1
               END IF
            END IF
 170     CONTINUE

      else

         DO 270 I2 = 1, IN2SV
            INP2 = IX2(I2)
            if (inp2 .le. 0) goto 270
            X = XN2(INP2)
            Y = YN2(INP2)
            IF (NDIM .GE. 3) Z = ZN2(INP2)
            DMIN = 1.0E38
            NDMIN = 0
            I1BEGS = I1BEG
            DO 260 I1 = I1BEGS, IN1SV
               NCOMP = NCOMP + 1
               nout  = nout  + 1
               IF ( nout .ge. 100000) THEN
                  NCMPX = (IN2SV - I2) * IN1 + (IN1 - I1)
                  WRITE (*, 10020) NCOMP, NCMPX
                  nout = 0
               END IF
               INP1 = IX1(I1)
               if (inp1 .le. 0) goto 260
               IF (NDIM .GE. 3) Z3D = ZN1(INP1)

               if ((imat .eq. 1 .and. (x-eps .gt. xn1(inp1))) .or.
     *           (imat .eq. 2 .and. (y-eps .gt. yn1(inp1))) .or.
     *           (imat .eq. 3 .and. (z-eps .gt. z3d))) then
                 i1beg = i1
                 goto 260
               end if

C     ... Since we are sorted on coordinat X|Y|Z,
C     if set 1 X|Y|Z greater than set 2 X|Y|Z+eps, go to next X2|Y2|Z2 coord.
               if ((imat .eq. 1 .and. xn1(inp1)-eps .gt. x) .or.
     *           (imat .eq. 2 .and. yn1(inp1)-eps .gt. y) .or.
     *           (imat .eq. 3 .and. z3d      -eps .gt. z)) goto 265

               DIS = MAX (ABS (XN1(INP1) - X), ABS(YN1(INP1) - Y),
     *              ABS (Z3D - Z) )
               IF ( (DIS .LE. EPS .AND. .NOT. CLOSE) .OR.
     $              DIS .EQ. 0.0) THEN
                  DISMAX = MAX(DISMAX, DIS)
                  NMATCH = NMATCH + 1
                  IXNP2(INP2) = INP1
                  IX2(I2) = -IX2(I2)
                  IN2 = IN2 - 1
                  IX1(I1) = 0
                  IN1 = IN1 - 1
                  GOTO 270
               ELSE IF (DIS .LE. EPS .AND. CLOSE) THEN
                  IF (DIS .LT. DMIN) THEN
                     DMIN = DIS
                     NDMIN = I1
                  END IF
               ELSE
                  DISMIN = MIN(DISMIN, DIS)
               END IF
 260        CONTINUE
 265        CONTINUE
            IF (CLOSE) THEN
               IF (DMIN .LE. EPS) THEN
                  DISMAX = MAX(DISMAX, DMIN)
                  NMATCH = NMATCH + 1
                  INP1 = IX1(NDMIN)
                  IXNP2(INP2) = INP1
                  IX2(I2) = -IX2(I2)
                  IN2 = IN2 - 1
                  IX1(NDMIN) = 0
                  IN1 = IN1 - 1
               END IF
            END IF
 270     CONTINUE

      end if

      CALL EXCPUS(time2)
      IF (BYSET) THEN
         IF ((IN1 .GT. 0) .AND. (IN2 .GT. 0)) THEN
            CALL PRTERR ('WARNING',
     &           'All nodes in nodal point set cannot be matched')
         END IF
      END IF

      WRITE (*, 10050) NCOMP
      WRITE (*, 10021) EPS
      IF (DISMAX .GT. -1.0E37) THEN
         WRITE (*, 10025) DISMAX
      END IF
      IF (DISMIN .LT. 1.0E37) THEN
         WRITE (*, 10030) DISMIN
      END IF
      IF (CLOSE) THEN
         WRITE (*, 10035)
      END IF
 180  CONTINUE
      IF (INIT) THEN
         WRITE (STRING, 10000) NMATCH
         CALL SQZSTR (STRING, LSTR)
         WRITE (*, 10010) STRING(:LSTR)
         IF (IEXPCT .NE. 0 .AND. IEXPCT .NE. NMATCH) then
            call prterr ('ERROR',
     $           'Expected number of matches not equal to actual')
            STOP 'Node Match Error'
         END IF
      ELSE
         WRITE (STRING, 10000) NMATCH, NMATCH-NSVMAT
10000    FORMAT (I8, ' nodes matched', :, ', ', I8, ' this set')
         CALL SQZSTR (STRING, LSTR)
         WRITE (*, 10010) STRING(:LSTR)
         IF (IEXPCT .NE. 0 .AND. IEXPCT .NE. NMATCH-NSVMAT) then
            call prterr ('ERROR',
     $           'Expected number of matches not equal to actual')
            STOP 'Node Match Error'
         END IF
      END IF
      if (ncomp .gt. 0) then
         if ((time2 - time1) .ne. 0.0) then
            Write (*, 10055) Time2-Time1, NCOMP/(time2-time1)
         else
            write (*, 10060) Time2-Time1
         end if
      end if

      RETURN
10010 FORMAT (/, 4X, A)
10020 FORMAT (' Number of equivalence comparisons = ',T38,I15,T60,I15)
10050 FORMAT (' Number of equivalence comparisons = ',T60,I15)
10021 FORMAT (' Tolerance used for matching = ',T60,1PE10.3)
10025 FORMAT (' Maximum distance between matched nodes    = ',T60,
     *     1PE10.3)
10030 FORMAT (' Minimum distance between nonmatched nodes = ',T60,
     *     1PE10.3)
10035 FORMAT (' Equivalencing based on closest node within tolerance')
10055 FORMAT (' Cpu Time = ',1PE10.3, ', comparison/sec = ',1PE10.3)
10060 FORMAT (' Cpu Time = ',1PE10.3, ', comparison/sec = Infinite')
      END
