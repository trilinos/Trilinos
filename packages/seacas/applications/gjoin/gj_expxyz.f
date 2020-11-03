C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C     -*- Mode: fortran -*-
C=======================================================================
      SUBROUTINE EXPXYZ (NDIM,
     $     MATNS1, MATNS2, NNNPS, IXNNPS, LTNNPS,
     $     NUMNP1, XN1, YN1, ZN1,
     $     NUMNP2, XN2, YN2, ZN2,
     $     IX1, IX2, nod1, nod2, IXNP2,
     $     nelbl1, ID1, nelb1, nlnk1, link1,
     $     nelbl2, id2, nelb2, nlnk2, link2,
     $     NMATCH, TOLER, CLOSE, MATMAT)
C=======================================================================

C     --*** EXPXYZ *** (GJOIN) Find matching nodes
C     --   Written by Greg Sjaardema, 8-25-92 - gdsjaar -
C     --   Modified from MATXYZ
C     --
C     --EXPXYZ matches nodes from the first and second database by comparing
C     --the coordinates.  Since some slack is needed in the equality check,
C     --the comparison may not work for all meshes.  The nodes to be matched
C     --may be limited to nodes in two nodal point sets. Nodes are limited
C     --by material blocks -- only match if the materials match.  Used to
C     --avoid matches across slidelines.
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
C     --   nod1 - SCRATCH - size = NUMNP1
C     --   nod2 - SCRATCH - size = NUMNP2
C     --   IXNP2 - IN/OUT - the index of the matching first database node;
C     --      negative if no match; reset if NMATCH = 0
C     --   NMATCH - IN/OUT - the number of matching nodes
C     --   CLOSE - IN/OUT - true if match closest node, false if match any
C     --      node within tolerance.
C     --   MATMAT - IN/OUT - true if match by material

      include 'gj_filnum.blk'
      INTEGER NNNPS(*), IXNNPS(*), LTNNPS(*)
      REAL XN1(*), YN1(*), ZN1(*)
      REAL XN2(*), YN2(*), ZN2(*)
      INTEGER IX1(*), IX2(*), nod1(*), nod2(*)
      INTEGER IXNP2(*)
      LOGICAL CLOSE, MATMAT

      INTEGER ID1(*), NELB1(*), NLNK1(*), LINK1(*)
      INTEGER ID2(*), NELB2(*), NLNK2(*), LINK2(*)

      PARAMETER (MXNAM = 1)
      DIMENSION KV(MXNAM), RV(MXNAM), IVAL(MXNAM)
      CHARACTER*8  CV(MXNAM)

      CHARACTER*80 STRING
      LOGICAL INIT
      LOGICAL BYSET
      LOGICAL OK

      NSVMAT = 0
      BYSET = (MATNS1 .GT. 0) .AND. (MATNS2 .GT. 0)

C     --Index nodes to match in the nodal point sets
      call prterr('WARNING',
     &     'EXPXYZ (By-material matching) is an experimental routine')
      call prterr('WARNING',
     &     'Use with caution and check validity of output and mesh')
      ioff1 = 0
      do 200 iblk = 1, nelbl1
         INIT = (NMATCH .LE. 0)
         id = id1(iblk)
         ioff2 = 0
         do 20 iblk2 = 1, nelbl2
            if (id .eq. id2(iblk2)) go to 30
            ioff2 = ioff2 + (nelb2(iblk2) * nlnk2(iblk2))
 20      continue
C     ... Fell through do loop - no matching material in mesh 2.
         write (*,*) 'Material ',id,' not present in mesh 2'
         go to 190
 30      continue

C     Now know that id1(iblk) == id2(iblk2).  Set flag in nod1 and nod2
C     corresponding to whether node exists in these blocks

         call iniint (numnp1, 0, nod1)
         call iniint (numnp2, 0, nod2)
         do 40 i = 1, (nelb1(iblk) * nlnk1(iblk))
            nod1(link1(i+ioff1)) = 1
 40      continue
         do 50 i = 1, (nelb2(iblk2) * nlnk2(iblk2))
            nod2(link2(i+ioff2)) = 1
 50      continue

         IF (BYSET) THEN
            IN1 = 0
            IX0 = IXNNPS(MATNS1) - 1
            DO 60 N = 1, NNNPS(MATNS1)
               if (nod1(ltnnps(ix0+n)) .eq. 1) then
                  IN1 = IN1 + 1
                  IX1(IN1) = LTNNPS(IX0+N)
               end if
 60         CONTINUE

            IN2 = 0
            IX0 = IXNNPS(MATNS2) - 1
            DO 70 N = 1, NNNPS(MATNS2)
               if (nod2(ltnnps(ix0+n)) .eq. 1) then
                  IN2 = IN2 + 1
                  IX2(IN2) = LTNNPS(IX0+N)
               end if
 70         CONTINUE
         else
C     ... Construct ix1 and ix2 from nod1 and nod2
            in1 = 0
            do 80 n = 1, numnp1
               if (nod1(n) .eq. 1) then
                  in1 = in1 + 1
                  ix1(in1) = n
               end if
 80         continue

            in2 = 0
            do 90 n = 1, numnp2
               if (nod2(n) .eq. 1) then
                  in2 = in2 + 1
                  ix2(in2) = n
               end if
 90         continue

         END IF

         time0 = 0.0
         time1 = 0.0
         time2 = 0.0

C     --Find the limits of the overlapping area of the two databases

         CALL MINMXS (IN1, IX1, XN1, X1MIN, X1MAX)
         CALL MINMXS (IN2, IX2, XN2, X2MIN, X2MAX)
         CALL MINMXS (IN1, IX1, YN1, Y1MIN, Y1MAX)
         CALL MINMXS (IN2, IX2, YN2, Y2MIN, Y2MAX)
         IF (NDIM .GE. 3) THEN
            CALL MINMXS (IN1, IX1, ZN1, Z1MIN, Z1MAX)
            CALL MINMXS (IN2, IX2, ZN2, Z2MIN, Z2MAX)
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

         EPS = ((XMAX - XMIN) + (YMAX - YMIN) + (ZMAX - ZMIN)) / 1E3
         IF (EPS .LT. 0.0) THEN
            CALL PRTERR ('ERROR', 'Nodes do not overlap')
            GOTO 180
         END IF

         IF (TOLER .GE. 0.0) THEN
            EPS = TOLER
         ELSE
            WRITE (*, 10040) EPS
10040       FORMAT (/' Default tolerance = ',1PE10.3)
            CALL FREFLD (0, 0,
     *           'Enter new value for tolerance (<ret> for default): ',
     *           MXNAM, IOS, NF, KV, CV, IVAL, RV)
            IF (IOS .NE. 0) RV(1) = 0.0
            IF (RV(1) .NE. 0.0) EPS = RV(1)
            CALL OUTLOG (KLOG, 1, 1, CV, IVAL, EPS)
         ENDIF
         OK = .TRUE.
         OK = OK .AND.
     &        ((MIN (X1MAX, X2MAX) - MAX (X1MIN, X2MIN)) .GT. -EPS)
         OK = OK .AND.
     &        ((MIN (Y1MAX, Y2MAX) - MAX (Y1MIN, Y2MIN)) .GT. -EPS)
         IF (NDIM .EQ. 3) THEN
            OK = OK .AND.
     &           ((MIN (Z1MAX, Z2MAX) - MAX (Z1MIN, Z2MIN)) .GT. -EPS)
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

         in1sv = in1
         IN1 = 0
         Z3D = 0.0
         DO 120 INP = 1, in1sv
           IF (IX1(INP) .GT. 0) THEN
             in = ix1(inp)
             IF ((XN1(IN) .GE. XMIN) .AND. (XN1(IN) .LE. XMAX)) THEN
               IF ((YN1(IN) .GE. YMIN) .AND. (YN1(IN) .LE. YMAX)) THEN
                 IF (NDIM .GE. 3) Z3D = ZN1(IN)
                 IF ((Z3D .GE. ZMIN) .AND. (Z3D .LE. ZMAX)) THEN
                   IN1 = IN1 + 1
                   IX1(IN1) = IN
                 END IF
               END IF
             END IF
           end if
 120       CONTINUE

           in2sv = in2
           IN2 = 0
           IF (NDIM .LT. 3) Z3D = 0.0
           DO 130 INP = 1, in2sv
             IF (IX2(INP) .GT. 0) THEN
               in = ix2(inp)
               IF ((XN2(IN) .GE. XMIN) .AND. (XN2(IN) .LE. XMAX)) THEN
                 IF ((YN2(IN) .GE. YMIN) .AND. (YN2(IN) .LE. YMAX)) THEN
                   IF (NDIM .GE. 3) Z3D = ZN2(IN)
                   IF ((Z3D .GE. ZMIN) .AND. (Z3D .LE. ZMAX)) THEN
                     IN2 = IN2 + 1
                     IX2(IN2) = IN
                   END IF
                 END IF
               END IF
             end if
 130     CONTINUE

C     --Blank out the IXNP2 array

         IF (INIT) THEN
            DO 140 INP = 1, NUMNP2
               IXNP2(INP) = -999
 140        CONTINUE
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
 150        CONTINUE
         END IF

C     --Find all matching nodes by comparing coordinates of nodes in overlap area

         call excpus(time0)
         call prterr('CMDSPEC', 'Entering Sorting Phase')
         call indexn(xn1, 1, 1, ix1, in1, .FALSE.)
         call indexn(xn2, 1, 1, ix2, in2, .FALSE.)
         call prterr('CMDSPEC', 'Entering Comparison Phase')

         time2 = 0.0
         CALL EXCPUS(time1)
         DISMIN =  1.0E38
         DISMAX = -1.0E38
         NCOMP  = 0
         nout   = 0
         WRITE (*,'(A)') ' '
         WRITE (*,'(A)') ' '

         NSVMAT = NMATCH
         Z3D = 0.0
         Z = 0.0
         IN1SV = IN1
         in2sv = in2
         i2beg = 1
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
               nout = nout + 1
               IF ( nout .ge. 100000) THEN
                  NCMPX = (IN1SV - I1) * IN2 + (IN2 - I2)
                  WRITE (*, 10020) NCOMP, NCMPX
                  nout = 0
               END IF
               INP2 = IX2(I2)
               if (inp2 .le. 0) goto 160
               if (x-eps .gt. xn2(inp2)) i2beg = i2
C     ... Since we are sorted on X, if set 2 X greater than set 1 X+eps,
C     go to next X1 coord.
               if (xn2(inp2)-eps .gt. x) goto 165

               IF (NDIM .GE. 3) Z3D = ZN2(INP2)
               DIS = MAX (ABS (XN2(INP2) - X), ABS(YN2(INP2) - Y),
     &              ABS (Z3D - Z) )
               IF ( (DIS .LE. EPS .AND. .NOT. CLOSE) .OR.
     &              DIS .EQ. 0.0) THEN
                  DISMAX = MAX(DISMAX, DIS)
                  NMATCH = NMATCH + 1
                  IXNP2(INP2) = INP1
                  IX2(I2) = -IX2(IN2)
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

         CALL EXCPUS(time2)
         IF (BYSET) THEN
            IF ((IN1 .GT. 0) .AND. (IN2 .GT. 0)) THEN
               CALL PRTERR ('WARNING',
     &              'All nodes in nodal point set cannot be matched')
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
 180     CONTINUE
         IF (INIT) THEN
            WRITE (STRING, 10000) ID, NMATCH
         ELSE
            WRITE (STRING, 10000) ID, NMATCH, NMATCH-NSVMAT
10000       FORMAT ('Material ID ',I8,' -->',
     &           I8, ' nodes matched', :, ', ', I8, ' this set')
         END IF
         CALL SQZSTR (STRING, LSTR)
         WRITE (*, 10010) STRING(:LSTR)
         if (ncomp .gt. 0) then
            if ((time2 - time1) .ne. 0.0) then
               Write (*, 10055) Time2-Time1, NCOMP/(time2-time1)
            else
               write (*, 10060) Time2-Time1
            end if
         end if
 190     continue
         ioff1 = ioff1 + (nelb1(iblk) * nlnk1(iblk))
 200  continue

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
