C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MATXYZ (NDIM,
     &     NUMNP, XN, YN, ZN, IX, IXNP, NMATCH, TOLER)
C=======================================================================

C     --*** MATXYZ *** (GREPOS) Find matching nodes
C     --
C     --MATXYZ matches nodes in the database by comparing
C     --the coordinates.  Since some slack is needed in the equality check,
C     --the comparison may not work for all meshes.
C     --The nodes to be matched may be limited to nodes in two
C     -- nodal point sets. (NOT IMPLEMENTED)
C     --
C     --Parameters:
C     --   NDIM - IN - the number of coordinates per node
C     --   NUMNP - IN - the number of nodes in the database
C     --   XN, YN, ZN - IN - the coordinates in the database
C     --   IX - SCRATCH - size = NUMNP
C     --   NMATCH - OUT - the number of matching nodes
C     --   TOLER - IN - the tolerance used for matching

      REAL XN(*), YN(*), ZN(*)
      INTEGER IX(*)
      INTEGER IXNP(*)

      CHARACTER*80 STRING

      IF (TOLER .LT. 0.0) THEN
         RETURN
      END IF

      nmatch = 0

      time0 = 0.0
      time1 = 0.0
      time2 = 0.0

      CALL MINMAX (NUMNP, XN, XMIN, XMAX)
      CALL MINMAX (NUMNP, YN, YMIN, YMAX)
      IF (NDIM .GE. 3) THEN
         CALL MINMAX (NUMNP, ZN, ZMIN, ZMAX)
      END IF

      DELTAX = XMAX - XMIN
      DELTAY = YMAX - YMIN
      DELTAZ = ZMAX - ZMIN

      XMIN = XMIN - TOLER
      XMAX = XMAX + TOLER
      YMIN = YMIN - TOLER
      YMAX = YMAX + TOLER
      ZMIN = ZMIN - TOLER
      ZMAX = ZMAX + TOLER

C     --Index the nodes

      IN1 = 0
      IF (NDIM .LT. 3) Z3D = 0.0
      DO 120 INP = 1, NUMNP
        IX(INP) = INP
        IXNP(INP) = INP
 120  CONTINUE

C     --Find all matching nodes by comparing coordinates of nodes in overlap area

      DELMAX = MAX(DELTAX, DELTAY, DELTAZ)

      imat = 0
      call excpus(time0)
      if (DELMAX .EQ. DELTAX) THEN
        imat = 1
        call indexn(xn, 1, 1, ix, numnp, .FALSE.)
        call prterr('CMDSPEC', 'Entering Sorting Phase, Sort on X')
      else if (delmax .eq. deltay) then
        imat = 2
        call indexn(yn, 1, 1, ix, numnp, .FALSE.)
        call prterr('CMDSPEC', 'Entering Sorting Phase, Sort on Y')
      else if (ndim .ge. 3 .and. delmax .eq. deltaz) then
        imat = 3
        call indexn(zn, 1, 1, ix, numnp, .FALSE.)
        call prterr('CMDSPEC', 'Entering Sorting Phase, Sort on Z')
      end if
      call prterr('CMDSPEC', 'Entering Comparison Phase')

      IN1 = NUMNP
      IN2 = NUMNP

      CALL EXCPUS(time1)
      DISMIN =  1.0E38
      DISMAX = -1.0E38
      NCOMP  = 0
      nout   = 0
      WRITE (*,'(A)') ' '
      WRITE (*,'(A)') ' '

      IF (NDIM .LT. 3) Z3D = 0.0
      IF (NDIM .LT. 3) Z = 0.0
      IN1SV = IN1
      IN2SV = IN2
      DO 170 I1 = 1, IN1SV
        I2BEG = I1+1
        INP1 = IX(I1)
        X = XN(INP1)
        Y = YN(INP1)
        IF (NDIM .GE. 3) Z = ZN(INP1)
        DMIN = 1.0E38
        NDMIN = 0

        I2BEGS = I2BEG
        DO 160 I2 = I2BEGS, IN2SV
          NCOMP = NCOMP + 1
          nout  = nout  + 1
          IF ( nout .ge. 10000000) THEN
             CMPX = 1.0 * (IN1SV - I1) * IN2 + (IN2 - I2)
             WRITE (*, 10020) NCOMP, CMPX, NMATCH
             nout = 0
          END IF
          INP2 = IX(I2)
          IF (NDIM .GE. 3) Z3D = ZN(INP2)

          if ((imat .eq. 1 .and. (x-toler .gt. xn(inp2))) .or.
     *        (imat .eq. 2 .and. (y-toler .gt. yn(inp2))) .or.
     *        (imat .eq. 3 .and. (z-toler .gt. z3d))) i2beg = i2

C ... Since we are sorted on coordinate X|Y|Z,
C     if set 2 X|Y|Z greater than set 1 X|Y|Z+toler, go to next X1|Y1|Z1 coord.
          if ((imat .eq. 1 .and. xn(inp2)-toler .gt. x) .or.
     *        (imat .eq. 2 .and. yn(inp2)-toler .gt. y) .or.
     *        (imat .eq. 3 .and. z3d     -toler .gt. z)) goto 165

          DIS = MAX (ABS (XN(INP2) - X), ABS(YN(INP2) - Y),
     *         ABS (Z3D - Z) )
          IF (DIS .EQ. 0.0) THEN
             DISMAX = MAX(DISMAX, DIS)
             NMATCH = NMATCH + 1
             if (inp1 .lt. inp2) then
                IXNP(INP2) = -INP1
             else
                IXNP(INP1) = -INP2
             end if
          ELSE IF (DIS .LE. TOLER) THEN
             IF (DIS .LT. DMIN) THEN
                DMIN = DIS
                NDMIN = I2
             END IF
          ELSE
             DISMIN = MIN(DISMIN, DIS)
          END IF
 160     CONTINUE
 165     CONTINUE
         IF (DMIN .LE. TOLER) THEN
            DISMAX = MAX(DISMAX, DMIN)
            NMATCH = NMATCH + 1
            INP2 = IX(NDMIN)
            if (inp1 .lt. inp2) then
               IXNP(INP2) = -INP1
            else
               IXNP(INP1) = -INP2
            end if
         END IF
 170   CONTINUE

      CALL EXCPUS(time2)

      WRITE (*, 10050) NCOMP, NMATCH
      WRITE (*, 10021) TOLER
      IF (DISMAX .GT. -1.0E37) THEN
         WRITE (*, 10025) DISMAX
      END IF
      IF (DISMIN .LT. 1.0E37) THEN
         WRITE (*, 10030) DISMIN
      END IF

      WRITE (STRING, 10000) NMATCH
      CALL SQZSTR (STRING, LSTR)
      WRITE (*, 10010) STRING(:LSTR)

      if (ncomp .gt. 0) then
         if ((time2 - time1) .ne. 0.0) then
            Write (*, 10055) Time2-Time1, NCOMP/(time2-time1)
         else
            write (*, 10060) Time2-Time1
         end if
      end if

      RETURN
10000 FORMAT (I8, ' nodes matched')
10010 FORMAT (/, 4X, A)
10020 FORMAT (' Comparisons = ',T20,I10,T32,1pE10.3,T47,
     *  ' Matches = ', I8)
10050 FORMAT (' Comparisons = ',T20,I10,T47,' Matches = ',I8)
10021 FORMAT (/,' Tolerance used for matching = ',T50,1PE10.3)
10025 FORMAT (' Maximum distance between matched nodes    = ',T50,
     *     1PE10.3)
10030 FORMAT (' Minimum distance between nonmatched nodes = ',T50,
     *     1PE10.3)
10055 FORMAT (' Cpu Time = ',1PE10.3, ', comparison/sec = ',1PE10.3)
10060 FORMAT (' Cpu Time = ',1PE10.3, ', comparison/sec = Infinite')
      END
