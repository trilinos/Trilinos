C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE HIDDEN (A, IHIDOP, LENF, NLNKE, NLNKF, LINKF,
     &   XN, YN, ZN, LENL, LINSET,
     &   MSCTYP, ZMMESH, MINMSH, MAXMSH, IELBST,
     &   DODEAD, IDN2B, HIDEF, HIDENP,
     &   ZF, DOIXF, NXFAC, IXFAC, NAMELB)
C=======================================================================

C   --*** HIDDEN *** (MESH) Hidden line removal routine
C   --
C   --HIDDEN determines if the face and nodes on a 3D surface
C   --are hidden.  The complexity of the algorithm used is determined by
C   --an option as follows:
C   --   1) A face is hidden if and only if its outward normal
C   --      points "into" the plotting surface.
C   --   2) A node is hidden if it is behind any surface.  A face
C   --      is hidden only if any of its nodes are hidden.
C   --   3) A node is hidden if it is behind any surface.  A face
C   --      is hidden only if any of its nodes are hidden.  An incomplete
C   --      algorithm creates partial line for mesh drawing.
C   --
C   --Dynamic memory is manipulated by the routine and should be checked
C   --for errors after call.
C   --
C   --Parameters:
C   --   A - IN/OUT - the dynamic memory base array
C   --   IHIDOP - IN - the hidden option
C   --   LENF - IN - the cumulative face counts by element block
C   --   NLNKE - IN - the number of nodes per element
C   --   NLNKF - IN - the number of nodes per face
C   --   LINKF - IN - the connectivity for all faces
C   --   XN, YN, ZN - IN/OUT - the nodal coordinates; the Z-coordinate is
C   --      pointing towards the viewer (out of the plotting plane)
C   --   LENL - IN - the cumulative line counts by element block
C   --   LINSET - IN/OUT - the sorted line set
C   --   MSCTYP - IN - the scale type that set ZMMESH (as in /MSHLIM/)
C   --   ZMMESH - IN - the zoom window
C   --   MINMSH, MAXMSH - IN - the minimum and maximum mesh lines needed
C   --      (as in MSHLIN of /MSHOPT/)
C   --   IELBST - IN - the element block status:
C   --      -1 = OFF, 0 = ON, but not selected, 1 = selected
C   --   DODEAD - IN - examine IDN2B iff true
C   --   IDN2B - IN - the element block for each dead node; dead if >= 0
C   --   HIDEF(i) - OUT - true iff face i is hidden
C   --   HIDENP(i) - OUT - true iff node i is hidden
C   --   ZF - IN - the face Z coordinates
C   --   DOIXF - IN - compute IXFAC iff true
C   --   NXFAC - OUT - the number of ordered faces (if DOIXF)
C   --   IXFAC - OUT - the indices of the ordered faces (if DOIXF)
C   --
C   --Common Variables:
C   --   Uses NELBLK of /DBNUMS/
C   --   Uses NUMNPF of /D3NUMS/

      PARAMETER (MSHNON=0, MSHBOR=1, MSHDIV=2, MSHSEL=3, MSHALL=4)
      PARAMETER (KLFT=1, KRGT=2, KBOT=3, KTOP=4, KNEA=5, KFAR=6)

      PARAMETER (KFVIS=0, KFNODH=10, KFPOUT=20, KFOUT=90, KFAWAY=100)
C      --Temporary face status:
C      --    0 = KFVIS  = visible
C      --   10 = KFNODH = some nodes hidden (set at end of processing)
C      --   20 = KFPOUT = outside zoom window, but within expanded limits
C      --   90 = KFOUT  = outside expanded limits
C      --  100 = KFAWAY = normal away
      PARAMETER (KNVIS=0, KNFOVR=10, KNHID=100)
C      --Temporary node status:
C      --   -n =          slide line with node n
C      --    0 = KNVIS  = visible
C      --   10 = KNFOVR = hidden by face
C      --  100 = KNHID  = hidden (normal away or outside zoom window)

      include 'params.blk'
      include 'debug.blk'
      include 'dbnums.blk'
      include 'd3nums.blk'

      DIMENSION A(*)
      INTEGER LENF(0:NELBLK)
      INTEGER NLNKE(NELBLK)
      INTEGER NLNKF(NELBLK)
      INTEGER LINKF(*)
      REAL XN(NUMNPF), YN(NUMNPF), ZN(NUMNPF)
      INTEGER LENL(-2:NELBLK), LINSET(LLNSET,*)
      CHARACTER*(*) MSCTYP
      REAL ZMMESH(KTOP)
      LOGICAL DODEAD
      INTEGER IELBST(NELBLK)
      INTEGER IDN2B(NUMNPF)
      INTEGER HIDEF(*)
      INTEGER HIDENP(NUMNPF)
      REAL ZF(*)
      LOGICAL DOIXF
      INTEGER IXFAC(*)
      CHARACTER*(MXSTLN) NAMELB(*)

      LOGICAL ZAWAY
      LOGICAL ZOOM
      LOGICAL SPEC, DOEDG
      REAL ZMLIM(KTOP)
      DO 110 IELB = 1, NELBLK
         DO 100 IFAC = LENF(IELB-1)+1, LENF(IELB)
            HIDEF(IFAC) = KFVIS
  100    CONTINUE
  110 CONTINUE

C   --Eliminate surfaces definitely outside zoom window and find smallest
C   --window limits

      IF ((MSCTYP .EQ. 'ZOOM') .OR.
     &   ((IHIDOP .GE. 2) .AND. (MSCTYP .NE. 'EACH'))) THEN
         CALL HIDZM (ZOOM,
     &      ZMMESH(KLFT), ZMMESH(KRGT), ZMMESH(KBOT), ZMMESH(KTOP),
     &      ZMLIM(KLFT), ZMLIM(KRGT), ZMLIM(KBOT), ZMLIM(KTOP),
     &      LENF, NLNKF, LINKF, XN, YN, HIDEF)
      ELSE
         ZOOM = .FALSE.
         IF ((IHIDOP .GE. 2) .AND. (MSCTYP .EQ. 'EACH')) THEN
            CALL CPYREA (KTOP, ZMMESH, ZMLIM)
         END IF
      END IF

C *** Simple hidden surface removal ***

      CALL INIINT (NUMNPF, KNHID, HIDENP)

C   --Mark node as visible if in visible face

      IF (IHIDOP .GE. 1) THEN

C      --Make shell faces face forward

         DO 130 IELB = 1, NELBLK
            IF (NLNKE(IELB) .EQ. 4) THEN
               IXL0 = IDBLNK (IELB, 0, LENF, NLNKF) - 1
               DO 120 IFAC = LENF(IELB-1)+1, LENF(IELB)
                  IF (HIDEF(IFAC) .LT. KFOUT) THEN
                     IF (ZAWAY (NLNKF(IELB), LINKF(IXL0+1),
     &                  XN, YN, ZN, HIDENP)) THEN
                        I = LINKF(IXL0+2)
                        LINKF(IXL0+2) = LINKF(IXL0+4)
                        LINKF(IXL0+4) = I
                     END IF
                  END IF
                  IXL0 = IXL0 + NLNKF(IELB)
  120          CONTINUE
            END IF
  130    CONTINUE

C      --Mark shell and truss nodes as visible if in visible face

         DO 160 IELB = 1, NELBLK
            IF (NLNKE(IELB) .LE. 4) THEN
               IXL0 = IDBLNK (IELB, 0, LENF, NLNKF) - 1
               DO 150 IFAC = LENF(IELB-1)+1, LENF(IELB)
                  IF (HIDEF(IFAC) .LT. KFOUT) THEN
                     DO 140 ILINK = 1, NLNKF(IELB)
                        HIDENP(LINKF(IXL0+ILINK)) = KNVIS
  140                CONTINUE
                  END IF
                  IXL0 = IXL0 + NLNKF(IELB)
  150          CONTINUE
            END IF
  160    CONTINUE
      END IF

      IF (IHIDOP .LE. 0) THEN

C      --Mark node as visible if in visible face

         DO 190 IELB = 1, NELBLK
            IXL0 = IDBLNK (IELB, 0, LENF, NLNKF) - 1
            DO 180 IFAC = LENF(IELB-1)+1, LENF(IELB)
               IF (HIDEF(IFAC) .LE. KFVIS) THEN
                  DO 170 ILINK = 1, NLNKF(IELB)
                     HIDENP(LINKF(IXL0+ILINK)) = KNVIS
  170             CONTINUE
               END IF
               IXL0 = IXL0 + NLNKF(IELB)
  180       CONTINUE
  190    CONTINUE

      ELSE IF (IHIDOP .EQ. 1) THEN

C      --Mark face as hidden if normal points away
C      --Mark node as visible if in visible face
         nhid = 0
         DO 210 IELB = 1, NELBLK
C            IF (NLNKE(IELB) .GT. 4) THEN
            IF (NLNKE(IELB) .GE. 4) THEN
               IXL = IDBLNK (IELB, 0, LENF, NLNKF)
               DO 200 IFAC = LENF(IELB-1)+1, LENF(IELB)
                  IF (HIDEF(IFAC) .LT. KFOUT) THEN
                     IF (ZAWAY (NLNKF(IELB), LINKF(IXL),
     &                  XN, YN, ZN, HIDENP)) THEN
                        HIDEF(IFAC) = KFAWAY
                        nhid = nhid + 1
                     END IF
                  END IF
                  IXL = IXL + NLNKF(IELB)
  200          CONTINUE
            END IF
  210    CONTINUE
         if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1))
     &      write (*, '(1x,a,i5)') 'faces with normal away =', nhid

      ELSE IF (IHIDOP .GE. 2) THEN

C      --Mark face as hidden if normals for all nodes point away
C      --Mark node as visible if normal points forward

         nhid = 0
         DO 230 IELB = 1, NELBLK
C            IF (NLNKE(IELB) .GT. 4) THEN
            IF (NLNKE(IELB) .Ge. 4) THEN
               IXL = IDBLNK (IELB, 0, LENF, NLNKF)
               DO 220 IFAC = LENF(IELB-1)+1, LENF(IELB)
                  IF (HIDEF(IFAC) .LT. KFOUT) THEN
                     IF (NZAWAY (NLNKF(IELB), LINKF(IXL),
     &                  XN, YN, ZN, HIDENP) .EQ. NLNKF(IELB)) THEN
                        HIDEF(IFAC) = KFAWAY
                        nhid = nhid + 1
                     END IF
                  END IF
                  IXL = IXL + NLNKF(IELB)
  220          CONTINUE
            END IF
  230    CONTINUE
         if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1))
     &      write (*, '(1x,a,i5)') 'faces with all normals away =', nhid
      END IF

C   --Mark node as visible if "dead" node, may be hidden later

      IF (DODEAD) THEN
         DO 240 INP = 1, NUMNPF
            IF (IDN2B(INP) .GE. 0) HIDENP(INP) = KNVIS
  240    CONTINUE
      END IF

C *** Hidden surface removal by hidden nodes ***

      IF (IHIDOP .GE. 2) THEN

C      --Assign each node to a block in the mesh

         rnumnpf = numnpf
         maxcr = sqrt(rnumnpf)
         MAXCR = max(20, maxcr)
         maxcr = min(5000, maxcr)
         CRDELT = MAXCR / MIN
     &      (ZMLIM(KRGT) - ZMLIM(KLFT), ZMLIM(KTOP) - ZMLIM(KBOT))
         CREPS = .001 / CRDELT / MAXCR
         COLMIN = ZMLIM(KLFT) - CREPS
         ROWMIN = ZMLIM(KBOT) - CREPS
         NUMCOL = INT (((ZMLIM(KRGT)+CREPS) - COLMIN) * CRDELT)
         NUMROW = INT (((ZMLIM(KTOP)+CREPS) - ROWMIN) * CRDELT)

         CALL MDRSRV ('IXNCRS', KIXCRS, (NUMCOL+1)*(NUMROW+1))
         CALL MDRSRV ('IXNCRE', KIXCRE, (NUMCOL+1)*(NUMROW+1))
         CALL MDRSRV ('NPCR', KNPCR, NUMNPF)
         CALL MDRSRV ('NPX', KNPX, NUMNPF)
         CALL MDRSRV ('NPY', KNPY, NUMNPF)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 360

         CALL INPCR (NUMNPF, HIDENP, XN, YN, ZN,
     &      COLMIN, ROWMIN, CRDELT, NUMCOL, NUMROW,
     &      A(KIXCRS), A(KIXCRE), A(KNPCR), LNPCR, A(KNPX), A(KNPY))

         CALL MDDEL ('NPX')
         CALL MDDEL ('NPY')
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 360

C      --Mark node as hidden if behind any face

         nhid = 0
         DO 260 IELB = 1, NELBLK
            IF (NLNKE(IELB) .GE. 4) THEN
               IXL = IDBLNK (IELB, 0, LENF, NLNKF)
               DO 250 IFAC = LENF(IELB-1)+1, LENF(IELB)
                  IF (HIDEF(IFAC) .LT. KFOUT) THEN
                     CALL HIDNOD (NLNKF(IELB), LINKF(IXL), XN, YN, ZN,
     &                  COLMIN, ROWMIN, CRDELT, CREPS, NUMCOL, NUMROW,
     &                  A(KIXCRS), A(KIXCRE), A(KNPCR), HIDENP, nhid)
                  END IF
                  IXL = IXL + NLNKF(IELB)
  250          CONTINUE
            END IF
  260    CONTINUE
         if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1)) then
            nslide = 0
            do 270 i = 1, numnpf
               if (hidenp(i) .eq. -1) then
CDEBUG            write (99, *) 'slide node =', i
                  nslide = nslide + 1
               end if
  270       continue
            write (*, '(1x,a,i5)') 'slide nodes =', nslide
         end if

         if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1))
     &      write (*, '(1x,a,i5,10x,a,i5)')
     &      'nodes visible =', lNPCR-nhid, 'nodes hidden =', nhid
         CALL MDDEL ('NPCR')
         CALL MDDEL ('IXNCRS')
         CALL MDDEL ('IXNCRE')
      END IF

C *** Identify edge of mesh with hidden faces ***

      NEDGES = 0
      IF (IHIDOP .GE. 1) THEN

C      --Find out if edge processing can be eliminated or cut down

         CALL CNTELB (IELBST, NELBLK, NUMON, NUMSEL)

         IF (IHIDOP .GE. 3) THEN
            DOEDG = .TRUE.
            SPEC = .FALSE.
         ELSE
            IF (MINMSH .GT. MSHSEL) THEN
               DOEDG = .FALSE.
            ELSE IF (MINMSH .EQ. MSHSEL) THEN
               DOEDG = (NUMSEL .LT. NUMON)
               SPEC = .TRUE.
            ELSE
               DOEDG = .TRUE.
               SPEC = .FALSE.
            END IF
         END IF

         IF (DOEDG) THEN

C         --Find the lines forming the exterior of the visible faces

            CALL MDRSRV ('IEDSET', KIESET, 3 * LENL(NELBLK))
            CALL MDRSRV ('NREF', KNREF, NUMNPF)
            CALL MDRSRV ('LREF', KLREF, NUMNPF)
            CALL MDRSRV ('MREF', KMREF, NUMNPF)
            CALL MDRSRV ('NBACK', KNBACK, NUMNPF)
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 360
            CALL HIDEDG (HIDEF, HIDENP, LENF, NLNKE, NLNKF, LINKF,
     &         SPEC, IELBST, A(KIESET), NEDGES,
     &         A(KNREF), A(KLREF), A(KMREF), A(KNBACK), NAMELB)

            CALL MDDEL ('NREF')
            CALL MDDEL ('LREF')
            CALL MDDEL ('MREF')
            CALL MDDEL ('NBACK')
            IF (NEDGES .GT. 0) THEN
               CALL MDLONG ('IEDSET', KIESET, 3 * NEDGES)
            ELSE
               CALL MDDEL ('IEDSET')
            END IF
            CALL MDSTAT (NERR, MEM)
            IF (NERR .GT. 0) GOTO 360
         END IF
      END IF

C *** Initialize hidden/partial line set flag ***

      DO 280 IL = 1, LENL(NELBLK)
         IF ((HIDENP(LINSET(1,IL)) .GT. KNVIS) .OR.
     &      (HIDENP(LINSET(2,IL)) .GT. KNVIS)) THEN
            LINSET(3,IL) = 0
         ELSE
            LINSET(3,IL) = 1
         END IF
  280 CONTINUE

C *** Hidden surface removal by partial lines ***

      IF (IHIDOP .GE. 3) THEN

C      --Find the partial line set and order so that hidden nodes are in
C      --LINSET(2,x)

         MAXLIN = LENL(NELBLK)
         CALL MDRSRV ('IPSET', KIPSET, MAXLIN)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 360
         SPEC = MAXMSH .LT. MSHSEL
         CALL HIDSET (LINSET, LENL, HIDENP, SPEC, IELBST,
     &      A(KIESET), NEDGES, A(KIPSET), NPART)

C      --Eliminate partial lines that are outside the zoom window

         IF (ZOOM) THEN
            CALL ZMSET
     &         (ZMMESH(KLFT), ZMMESH(KRGT), ZMMESH(KBOT), ZMMESH(KTOP),
     &         XN, YN, LINSET, A(KIPSET), NPART)
         END IF

         CALL MDLONG ('IPSET', KIPSET, NPART)
         CALL MDRSRV ('LREF', KLREF, NUMNPF)
         CALL MDRSRV ('MREF', KMREF, NUMNPF)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 360

C      --Find line to corner point which is totally covered

         CALL HIDCOR (LINSET, A(KIPSET), NPART,
     &      A(KIESET), NEDGES, LENF, NLNKF, LINKF, XN, YN, ZN,
     *     A(KMREF), A(KLREF))

         CALL MDDEL ('MREF')
         CALL MDDEL ('LREF')

C      --Find crossing point for each partial line

         CALL MDRSRV ('TVHM', KTVHM, NPART)
         CALL MDRSRV ('ICROSS', KCROSS, NPART)
         CALL MDRSRV ('XP', KXP, NPART)
         CALL MDRSRV ('YP', KYP, NPART)
         CALL MDRSRV ('ZP', KZP, NPART)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 360
         CALL HIDPAR (LINSET, A(KIPSET), NPART,
     &      A(KIESET), NEDGES, LENF, NLNKE, NLNKF, LINKF,
     &      XN, YN, ZN, HIDEF, HIDENP,
     &      A(KTVHM), A(KCROSS), A(KXP), A(KYP), A(KZP))

         CALL MDDEL ('TVHM')
         CALL MDDEL ('ICROSS')
         CALL MDDEL ('XP')
         CALL MDDEL ('YP')
         CALL MDDEL ('ZP')
         CALL MDDEL ('IPSET')
      END IF

      IF (IHIDOP .GE. 3) THEN
         IF (NEDGES .GT. 0) THEN

C         --Find out if some edges can be eliminated

            IF (MINMSH .EQ. MSHSEL) THEN
               IF (NUMSEL .EQ. NUMON) THEN
                  NEDGES = 0
               ELSE
                  CALL MATEDG (LENF, IELBST, A(KIESET), NEDGES)
               END IF
               IF (NEDGES .LE. 0) THEN
                  CALL MDDEL ('IEDSET')
               END IF
            END IF
         END IF
      END IF

      IF (NEDGES .GT. 0) THEN

C      --Eliminate edges that are outside the zoom window

         IF (ZOOM) THEN
            CALL ZMEDGE (
     &         ZMMESH(KLFT), ZMMESH(KRGT), ZMMESH(KBOT), ZMMESH(KTOP),
     &         XN, YN, A(KIESET), NEDGES)
         END IF

C      --Identify edges within LINSET

         CALL MDRSRV ('IXLIN', KIXLIN, LENL(NELBLK))
         CALL MDRSRV ('MREF', KMREF, NUMNPF)
         CALL MDRSRV ('LREF', KLREF, NUMNPF)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 360

         CALL IDEDGE (A(KIESET), NEDGES, HIDENP,
     &      LINSET, LENL(NELBLK), A(KIXLIN),
     *      A(KMREF), A(KLREF))

         CALL MDDEL ('LREF')
         CALL MDDEL ('MREF')
         CALL MDDEL ('IXLIN')
         CALL MDDEL ('IEDSET')
      END IF

C   --Make up indexed list of faces (back to front)

      IF (DOIXF) THEN
         CALL MDRSRV ('SCRIX', KSCRIX, LENF(NELBLK))
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 360

         CALL HIDIXF (HIDEF, HIDENP, LENF, NLNKE, NLNKF, LINKF,
     &      ZF, NXFAC, IXFAC, A(KSCRIX), NAMELB)

         CALL MDDEL ('SCRIX')
      END IF

      IF (IHIDOP .GE. 2) THEN

C      --Mark face as hidden if any of its nodes are hidden

         nhid = 0
         DO 320 IELB = 1, NELBLK
            IXL0 = IDBLNK (IELB, 0, LENF, NLNKF) - 1
            DO 310 IFAC = LENF(IELB-1)+1, LENF(IELB)
               IF (HIDEF(IFAC) .LE. KFVIS) THEN
                  DO 290 ILINK = 1, NLNKF(IELB)
                     IF (HIDENP(LINKF(IXL0+ILINK)) .GT. KNVIS) THEN
                        HIDEF(IFAC) = KFNODH
                        nhid = nhid + 1
                        GOTO 300
                     END IF
  290             CONTINUE
  300             CONTINUE
               END IF
               IXL0 = IXL0 + NLNKF(IELB)
  310       CONTINUE
  320    CONTINUE
         if ((cdebug .eq. 'HIDDEN') .and. (idebug .ge. 1))
     &      write (*, '(1x,a,i5)') 'faces with hidden node =', nhid
      END IF

C   --Change the node and face status to logicals

      DO 330 INP = 1, NUMNPF
         IF (HIDENP(INP) .LE. KNVIS) THEN
            HIDENP(INP) = 0
         ELSE
            HIDENP(INP) = -1
         END IF
  330 CONTINUE

      DO 350 IELB = 1, NELBLK
         DO 340 IFAC = LENF(IELB-1)+1, LENF(IELB)
            IF (HIDEF(IFAC) .LE. KFVIS) THEN
               HIDEF(IFAC) = 0
            ELSE
               HIDEF(IFAC) = -1
            END IF
  340    CONTINUE
  350 CONTINUE

  360 CONTINUE
      RETURN
      END
