C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE HIDPAR (LINSET, IPSET, NPART, IEDSET, NEDGES,
     &           LENF, NLNKE, NLNKF, LINKF, XN, YN, ZN, HIDEF,
     &           HIDENP, TVHMAX, ICROSS, XP, YP, ZP)
C=======================================================================

C   --*** HIDPAR *** (MESH) Create new nodes for partial lines
C   --   Written by Amy Gilkey - revised 02/24/88
C   --
C   --HIDPAR finds the amount of a partial line that is hidden.  It then
C   --creates a new node at the last visible point and changes the LINSET(3,x)
C   --to point to the new node.
C   --
C   --Parameters:
C   --   LINSET   - I/O - the sorted line set
C   --   IPSET    - I/O - the indices of the partial line set
C   --   NPART    - I/O - the number of lines in the partial line set
C   --   IEDSET   - IN  - the edge line set;
C   --                    (0) = face defining edge; 0 to delete edge
C   --   NEDGES   - I/O - the number of lines in the edge set
C   --   LENF     - IN  - the cumulative face counts by element block
C   --   NLNKE    - IN  - the number of nodes per element
C   --   NLNKF    - IN  - the number of nodes per face
C   --   LINKF    - IN  - the connectivity for all faces
C   --   XN,YN,ZN - IN  - the nodal coordinates
C   --   HIDEF(i) - IN  - true iff face i is hidden
C   --   HIDENP   - I/O - node status (as in HIDDEN)
C   --   TVHMAX   - SCR - size is NPART
C   --   CROSS    - SCR - size is NPART
C   --
C   --Common Variables:
C   --   Uses NUMNPF, LLNSET of /D3NUMS/

      PARAMETER (KFVIS=0, KFNODH=10, KFPOUT=20, KFOUT=90, KFAWAY=100)
      PARAMETER (KNVIS=0, KNFOVR=10, KNHID=100)

      PARAMETER (EPS = .001)
C   --EPS is a normalized error allowance

      include 'debug.blk'
      include 'dbnums.blk'
      include 'd3nums.blk'
      include 'sizes.blk'

      INTEGER LINSET(LLNSET,*)
      INTEGER IPSET(*)
      INTEGER IEDSET(0:2,*)
      INTEGER LENF(0:NELBLK)
      INTEGER NLNKE(NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      REAL XN(*), YN(*), ZN(*)
      INTEGER HIDEF(*)
      INTEGER HIDENP(*)
      REAL TVHMAX(*)
      INTEGER ICROSS(*)
      REAL XP(*), YP(*), ZP(*)

      LOGICAL CKCROS

C   --Clear out partial line parameter

      CALL INIREA (NPART, -99.0, TVHMAX)
      CALL INIINT (NPART, 0, ICROSS)

C   --Find an epsilon based on the average edge size

      XLEN = 0.0
      YLEN = 0.0
      DO 100 IEDG = 1, NEDGES
         IF (IEDSET(0,IEDG) .EQ. 0) GOTO 100
         N1 = IEDSET(1,IEDG)
         N2 = IEDSET(2,IEDG)
         X1 = XN(N1)
         X2 = XN(N2)
         Y1 = YN(N1)
         Y2 = YN(N2)
         XLEN = XLEN + ABS (X2-X1)
         YLEN = YLEN + ABS (Y2-Y1)
  100 CONTINUE

      IF (NEDGES .GT. 0) THEN
         XLEN = XLEN / NEDGES
         YLEN = YLEN / NEDGES
         EPSDAT = SQRT (XLEN**2 + YLEN**2) * .01
      END IF

C   --Check each edge line against each partial line for overlap

      NOLDPT = NPART

      nhid = 0
      DO 130 IEDG = 1, NEDGES
         IFAC = IEDSET(0,IEDG)
         IF (IFAC .EQ. 0) GOTO 130
         N1 = IEDSET(1,IEDG)
         N2 = IEDSET(2,IEDG)

C      --Calculate X-Y-Z box enclosing edge line

         X1 = XN(N1)
         X2 = XN(N2)
         XMIN = MIN (X1, X2) - EPSDAT
         XMAX = MAX (X1, X2) + EPSDAT
         Y1 = YN(N1)
         Y2 = YN(N2)
         YMIN = MIN (Y1, Y2) - EPSDAT
         YMAX = MAX (Y1, Y2) + EPSDAT
         Z1 = ZN(N1)
         Z2 = ZN(N2)
         ZMIN = MIN (Z1, Z2) - EPSDAT
         ZMAX = MAX (Z1, Z2) + EPSDAT

         XLN = X2 - X1
         YLN = Y2 - Y1

         IP = 1
  110    CONTINUE
         IF (IP .LE. NPART) THEN
            IH = LINSET(2,IPSET(IP))
            IV = LINSET(1,IPSET(IP))

C         --Determine if partial line is within X-Y-Z box enclosing edge line

            X0 = XN(IH)
            XV = XN(IV)
            IF (XMAX .LT. MIN (X0, XV)) GOTO 120
            IF (XMIN .GT. MAX (X0, XV)) GOTO 120
            Y0 = YN(IH)
            YV = YN(IV)
            IF (YMAX .LT. MIN (Y0, YV)) GOTO 120
            IF (YMIN .GT. MAX (Y0, YV)) GOTO 120
            Z0 = ZN(IH)
            ZV = ZN(IV)
            IF (ZMAX .LT. MIN (Z0, ZV)) GOTO 120

            IF ((N1 .EQ. IV) .OR. (N2 .EQ. IV)) THEN
               IF (N2 .NE. IH) ICROSS(IP) = IEDG
               GOTO 120
            END IF
            IF (N2 .EQ. IH) GOTO 120

C         --Calculate the intersection of the edge and the partial line
C         --Solve the simultaneous equations:
C         --   X = X0 + (XV - X0) * TVH = X1 + (X2 - X1) * TLN
C         --   Y = Y0 + (YV - Y0) * TVH = Y1 + (Y2 - Y1) * TLN

            XVH = XV - X0
            YVH = YV - Y0
            XLH = X1 - X0
            YLH = Y1 - Y0
            DET = XVH * (-YLN) - YVH * (-XLN)
            IF (DET .EQ. 0.0) GOTO 120

            TVH = (-YLN * XLH + XLN * YLH) / DET
            IF ((TVH .GE. 0) .AND. (TVH .LE. 1+EPS)) THEN
               TLN = (-YVH * XLH + XVH * YLH) / DET
               IF ((TLN .GE. -EPS) .AND. (TLN .LE. 1+EPS)) THEN

C               --Save the overlap farthest from the hidden node

                  IF (TVHMAX(IP) .LT. TVH) THEN
                     IF (ZMIN .GT. Z0) THEN
                        TVHNEW = TVH
                     ELSE
                        ZCR = Z0 + (ZV - Z0) * TVH
                        ZCLN = Z1 + (Z2 - Z1) * TLN
                        IF (ZCR .LE. ZCLN) THEN
                           TVHNEW = TVH
                        ELSE
                           TVHNEW = TVHMAX(IP)
                        END IF
                     END IF

C                  --Delete totally hidden lines from partial set

                     IF (TVHNEW .GE. 1-EPS) THEN
                        IF (TLN .LE. 0.5) THEN
                           N = N1
                        ELSE
                           N = N2
                        END IF
                        IELB = 0
                        IXL = IDBLNK (IELB, IFAC, LENF, NLNKF)
                        IF (CKCROS (N, IH, NLNKF(IELB), LINKF(IXL),
     &                      XN, YN, ZN)) THEN

C                        --If face hides the line, move partial line to totally
C                        --hidden list

                           I = IPSET(IP)
                           IPSET(IP) = IPSET(NPART)
                           IPSET(NPART) = I
                           TVHMAX(IP) = TVHMAX(NPART)
                           ICROSS(IP) = ICROSS(NPART)
                           NPART = NPART - 1
                           IP = IP - 1
                        ELSE IF (ICROSS(IP) .EQ. 0) THEN
                           ICROSS(IP) = IEDG
                        END IF

                     ELSE IF ((TLN .GE. 0.0) .AND. (TLN .LE. 1.0)) THEN
                        TVHMAX(IP) = TVHNEW
                     END IF
                  END IF
               END IF
            END IF
  120       CONTINUE
            IP = IP + 1
            GOTO 110
         END IF
  130 CONTINUE
      if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1))
     &   write (*, '(1x,a,i5)') 'invisible lines =', noldpt-npart

C   --Group those partial lines which may end on an edge

      NNPART = NPART
      DO 140 IP = NNPART, 1, -1
         IF ((ICROSS(IP) .NE. 0) .OR. (TVHMAX(IP) .LT. EPS)) THEN
            I = IPSET(IP)
            IPSET(IP) = IPSET(NPART)
            IPSET(NPART) = I
            T = TVHMAX(IP)
            TVHMAX(IP) = TVHMAX(NPART)
            TVHMAX(NPART) = T
            I = ICROSS(IP)
            ICROSS(IP) = ICROSS(NPART)
            ICROSS(NPART) = I
            NPART = NPART - 1
         END IF
  140 CONTINUE
      if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1))
     &   write (*, '(1x,a,i5)') 'questionable lines =', NNPART-npart

C   --Find the midpoint of the questionable partial lines

      DO 150 IP = NPART+1, NNPART
         IH = LINSET(2,IPSET(IP))
         IV = LINSET(1,IPSET(IP))
         TVH = TVHMAX(IP)
         if (tvh .lt. 0) tvh = 0.0
         XP(IP) = 0.5 * (XN(IV) + XN(IH) + (XN(IV) - XN(IH)) * TVH)
         YP(IP) = 0.5 * (YN(IV) + YN(IH) + (YN(IV) - YN(IH)) * TVH)
         ZP(IP) = 0.5 * (ZN(IV) + ZN(IH) + (ZN(IV) - ZN(IH)) * TVH)
  150 CONTINUE

C   --Find out if the midpoint of the questionable partial lines are
C   --hidden by a visible face

      NQUES = NNPART - NPART
      IQUES = NPART+1
      DO 170 IELB = 1, NELBLK
         IF (NLNKE(IELB) .GE. 4) THEN
            IXL = IDBLNK (IELB, 0, LENF, NLNKF)
            DO 160 IFAC = LENF(IELB-1)+1, LENF(IELB)
               IF (HIDEF(IFAC) .LT. KFOUT) THEN
                  CALL HIDMID (NLNKF(IELB), LINKF(IXL),
     &               XN, YN, ZN, LINSET,
     &               IPSET(IQUES), TVHMAX(IQUES), ICROSS(IQUES),
     &               XP(IQUES), YP(IQUES), ZP(IQUES), NQUES)
               END IF
               IXL = IXL + NLNKF(IELB)
  160       CONTINUE
         END IF
  170 CONTINUE
      NPART = NPART + NQUES
      if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1))
     &   write (*, '(1x,a,i5)') 'invisible lines =', NNPART-npart

C   --Delete the edges which are totally hidden lines

      nhid = 0
      DO 190 IP = NPART+1, NOLDPT
         IV = LINSET(1,IPSET(IP))
         IH = LINSET(2,IPSET(IP))
         DO 180 IEDG = 1, NEDGES
            IF (IEDSET(1,IEDG) .EQ. IV) THEN
               IF (IEDSET(2,IEDG) .EQ. IH) THEN
                  IEDSET(0,IEDG) = 0
                  nhid = nhid + 1
                  GOTO 190
               END IF
            END IF
  180    CONTINUE
  190 CONTINUE
      if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1))
     &   write (*, '(1x,a,i5)') 'invisible edges =', nhid

C   --Delete totally visible lines from partial set

      nwhole = 0
      DO 200 IP = IQUES, NPART
         IF (TVHMAX(IP) .LT. 0.001) THEN
C         --Line is totally visible, so delete from partial line set as visible
C         --and mark hidden node as visible
            nwhole = nwhole + 1
            IH = LINSET(2,IPSET(IP))
            HIDENP(IH) = KNVIS
            LINSET(3,IPSET(IP)) = 1
         END IF
  200 CONTINUE
      if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1))
     &   write (*, '(1x,a,i5)') 'whole lines =', nwhole

C   --Calculate the visible part of each partial line, put coordinates into
C   --XP, YP, ZP (cannot overwrite coordinates until done using them
C   --for partial line calculation)

      DO 210 IP = 1, NPART
         IF (LINSET(3,IPSET(IP)) .EQ. 0) THEN
            IH = LINSET(2,IPSET(IP))
            IV = LINSET(1,IPSET(IP))
            TVH = TVHMAX(IP)
            XP(IP) = XN(IH) + (XN(IV) - XN(IH)) * TVH
            YP(IP) = YN(IH) + (YN(IV) - YN(IH)) * TVH
            ZP(IP) = ZN(IH) + (ZN(IV) - ZN(IH)) * TVH
         END IF
  210 CONTINUE

C   --Make new partial line nodes by filling in coordinates of hidden nodes
C   --and pointing LINSET(3,i) to these nodes

C   --Skip node 1 since LINSET(3,x) = 1 is reserved for whole lines
      IPART = NUMNPF

      DO 230 IP = 1, NPART
         IF (LINSET(3,IPSET(IP)) .EQ. 0) THEN
            IPART = IPART + 1
            IF (IPART .GT. NPSIZ) THEN
               WRITE (*, 10000) 'ERROR in HIDPAR', IP, IPART, NUMNPF
10000           FORMAT (1X, A, 3I5)
               GOTO 240
            END IF
            LINSET(3,IPSET(IP)) = IPART
            XN(IPART) = XP(IP)
            YN(IPART) = YP(IP)
            ZN(IPART) = ZP(IP)
         END IF
  230 CONTINUE
  240 CONTINUE
      if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1))
     &   write (*, '(1x,a,i5)') 'partial lines =', npart

C   --Reset partial line set to include totally KNVISible lines
      NPART = NNPART

      RETURN
      END
