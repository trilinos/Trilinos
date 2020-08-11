C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE QMESH (A, IA, MP, ML, MS, MR, MSC, MCOM, ICOM, JCOM,
     &   CIN, RIN, IIN, KIN, IUNIT, IDUMP, N, IPOINT, COOR, IPBOUN,
     &   ILINE, LTYPE, NINT, FACTOR, LCON, ILBOUN, ISBOUN, ISIDE, NLPS,
     &   IFLINE, ILLIST, IBARST, JMAT, JCENT, NLPB, JFLINE, JLLIST,
     &   IREGN, IMAT, NSPR, IFSIDE, ISLIST, IRPB, IPBF, NPPF, IFPB,
     &   LISTPB, ILBF, NLPF, IFLB, LISTLB, ISBF, NSPF, IFSB, LISTSB,
     &   LINKP, LINKL, LINKS, LINKB, LINKR, LINKSC, LINKPB, LINKLB,
     &   LINKSB, RSIZE, IFHOLE, NHPR, IHLIST, IRGFLG, ISCHM, SCHEME,
     &   DEFSCH, DEFSIZ, NPREGN, NPNBC, NPSBC, NPNODE, NPELEM, MAXKXN,
     &   STEP, DEV1, THREE, EIGHT, NINE, LGROUP, BATCH, AMESUR, XNOLD,
     &   YNOLD, NXKOLD, MMPOLD, LINKEG, LISTEG, BMESUR, MLINK, NPROLD,
     &   NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX,
     &   IDIVIS, SIZMIN, EMAX, EMIN)

C***********************************************************************

C  QMESH: A QUADRILATERAL MESH GENERATION PROGRAM

C***********************************************************************

C  ORIGINALLY WRITTEN BY:
C     RONDALL E JONES  DIV 2642  SANDIA LABORATORIES  ALBUQUERQUE
C  REWRITTEN AND UPDATED BY:
C     TEDDY D. BLACKER  DIV 1522 SANDIA LABORATORIES  ALBUQUERQUE
C     DECEMBER 1985

C***********************************************************************

      DIMENSION A(1), IA(1)
      DIMENSION IPOINT(MP), COOR(2, MP), IPBOUN(MP)
      DIMENSION ILINE(ML), LTYPE(ML), NINT(ML), FACTOR(ML), LCON(3, ML)
      DIMENSION ILBOUN(ML), ISBOUN(ML)
      DIMENSION ISIDE(MS), NLPS(MS), IFLINE(MS), ILLIST(MS*3)
      DIMENSION IBARST(MS), JMAT(MS), JCENT(MS), NLPB(MS), JFLINE(MS)
      DIMENSION JLLIST(MS*3)
      DIMENSION IREGN(MR), IMAT(MR), NSPR(MR), IFSIDE(MR), ISLIST(MR*4)
      DIMENSION IRPB(MR), RSIZE(MR), IFHOLE(MR), NHPR(MR), IHLIST(MR*2)
      DIMENSION ISCHM(MSC), SCHEME(MSC)
      DIMENSION IPBF(MP), NPPF(MP), IFPB(MP), LISTPB(2, MP)
      DIMENSION ILBF(ML), NLPF(ML), IFLB(ML), LISTLB(2, ML)
      DIMENSION ISBF(ML), NSPF(ML), IFSB(ML), LISTSB(2, ML)
      DIMENSION LINKP(2, MP), LINKL(2, ML), LINKS(2, MS), LINKB(2, MS)
      DIMENSION LINKR(2, MR), LINKSC(2, MR), LINKPB(2, MP)
      DIMENSION LINKLB(2, ML), LINKSB(2, ML), IRGFLG(MR)

      DIMENSION AMESUR(NPEOLD), XNOLD(NPNOLD), YNOLD(NPNOLD)
      DIMENSION NXKOLD(NNXK, NPEOLD), MMPOLD(3, NPROLD)
      DIMENSION LINKEG(2, MLINK), LISTEG(4 * NPEOLD), BMESUR(NPNOLD)

      DIMENSION N(29), CIN(MCOM), IIN(MCOM), RIN(MCOM), KIN(MCOM)

      DIMENSION K(30), IDUMMY(1)

      CHARACTER*72 SCHEME, DEFSCH, SCHSTR, CIN
      CHARACTER DEV1*3

      LOGICAL NOROOM, EVEN, ERR, CCW, IANS, LGROUP
      LOGICAL RECT, REAL, STEP, TEST, REMESH
      LOGICAL BAR, ADDLNK, EIGHT, NINE, PENTAG, TRIANG, TRNSIT, FINAL
      LOGICAL HALFC, COUNT, FILL, ERRCHK, THREE, BATCH, GRAPH

      DATA IEXIT, IOVER, IQUIT /1, 2, 3/

C  INITIALIZE

      IZ = 0
      IPNTR  = 0
      IPNTR1 = 0
      IPNTR2 = 0
      NPREGN = 0
      NPELEM = 1
      NPNODE = 1
      NPNBC = 1
      NPSBC = 1
      MAXKXN = 1
      MXND = 1
      MXNL = 1
      MXRXG = 1
      MXNPER = 1
      MAXSBC = 1
      MAXNBC = 1
      KSBC = 1
      KKSBC = 1
      MAX1 = 1
      MAX2 = 1
      MAX3 = 1
      MAX4 = 1
      EVEN = .TRUE.
      COUNT = .TRUE.
      ADDLNK = .FALSE.
      PENTAG = .FALSE.
      TRIANG = .FALSE.
      TRNSIT = .FALSE.
      FILL = .FALSE.
      GRAPH = .FALSE.

C  HEADER

      CALL MESAGE (' ')
      CALL MESAGE ('MESH PROCESSING BEGUN')
      CALL MESAGE (' ')

C  FILL IN ANY MISSING INTERVALS ACCORDING TO SIZE AND CHECK THE
C  VALIDITY OF REGION DATA

      ERRCHK = .FALSE.
      DO 130 I = 1, N(9)
         CALL LTSORT (MR, LINKR, ABS(IRPB(I)), IPNTR1, ADDLNK)
         CALL LTSORT (MS, LINKB, ABS(IRPB(I)), IPNTR2, ADDLNK)
         IF ((IRPB(I) .GT. 0) .AND. (IRPB(I) .LE. N(22)) .AND.
     &      (IPNTR1 .GT. 0)) THEN
            IF (IRGFLG(IPNTR1) .LE. -1) THEN
               L = IPNTR1
               IF (RSIZE (L) .LE. 0.) THEN
                  SIZE = DEFSIZ
               ELSE
                  SIZE = RSIZE (L)
               ENDIF
               CALL DATAOK (MP, ML, MS, MR, L, IREGN(L), COOR, ILINE,
     &            LTYPE, NINT, LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE,
     &            ISLIST, LINKP, LINKL, LINKS, SIZE, ERRCHK, ERR)
               IF (NHPR (L) .GT. 0) THEN
                  DO 100 IHOLE = IFHOLE(L), IFHOLE(L) + NHPR(L) - 1
                     CALL LTSORT (MR, LINKR, ABS(IHLIST(IHOLE)), IPNTRH,
     &                  ADDLNK)
                     IF (IPNTRH .GT. 0) THEN
                        LL = IPNTRH
                        IF (RSIZE (LL) .LE. 0.) THEN
                           SIZE = DEFSIZ
                        ELSE
                           SIZE = RSIZE (LL)
                        ENDIF
                        CALL DATAOK (MP, ML, MS, MR, LL, IREGN(LL),
     &                     COOR, ILINE, LTYPE, NINT, LCON, NLPS, IFLINE,
     &                     ILLIST, NSPR, IFSIDE, ISLIST, LINKP, LINKL,
     &                     LINKS, SIZE, ERRCHK, ERR)
                     ENDIF
  100             CONTINUE
               ENDIF
            ELSE IF (IRGFLG(IPNTR1) .GE. 1) THEN
               J1 = IFSIDE(IPNTR1)
               J2 = J1 + NSPR(IPNTR1) - 1
               IREGN(IPNTR1) = -ABS(IREGN(IPNTR1))
               MXRXG = MAX(MXRXG, NSPR(IPNTR1))
               DO 120 J = J1, J2
                  CALL LTSORT (MR, LINKR, ABS(ISLIST(J)), IPNTR1,
     &               ADDLNK)
                  IF (IPNTR1 .GT. 0) THEN
                     L = IPNTR1
                     IF (RSIZE (L) .LE. 0.) THEN
                        SIZE = DEFSIZ
                     ELSE
                        SIZE = RSIZE (L)
                     ENDIF
                     CALL DATAOK (MP, ML, MS, MR, L, IREGN(L), COOR,
     &                  ILINE, LTYPE, NINT, LCON, NLPS, IFLINE, ILLIST,
     &                  NSPR, IFSIDE, ISLIST, LINKP, LINKL, LINKS,
     &                  SIZE, ERRCHK, ERR)
                     IF (NHPR (L) .GT. 0) THEN
                        DO 110 IHOLE = IFHOLE(L),
     &                     IFHOLE(L) + NHPR(L) - 1
                           CALL LTSORT (MR, LINKR, ABS(IHLIST(IHOLE)),
     &                        IPNTRH, ADDLNK)
                           IF (IPNTRH .GT. 0) THEN
                              LL = IPNTRH
                              IF (RSIZE (LL) .LE. 0.) THEN
                                 SIZE = DEFSIZ
                              ELSE
                                 SIZE = RSIZE (LL)
                              ENDIF
                              CALL DATAOK (MP, ML, MS, MR, LL,
     &                           IREGN(LL), COOR, ILINE, LTYPE, NINT,
     &                           LCON, NLPS, IFLINE, ILLIST, NSPR,
     &                           IFSIDE, ISLIST, LINKP, LINKL, LINKS,
     &                           SIZE, ERRCHK, ERR)
                           ENDIF
  110                   CONTINUE
                     ENDIF
                  END IF
  120          CONTINUE
            END IF
         END IF
  130 CONTINUE
      ERRCHK = .TRUE.

C  FIND THE MAXIMUM NUMBER OF LINES/REGION, PERIMETER POINTS/REGION,
C  AND HOLES/REGION

      DO 140 I = 1, N(2)
         MAX1 = MAX0(NINT(I), MAX1)
  140 CONTINUE
      DO 150 I = 1, N(3)
         MAX2 = MAX0(NLPS(I), MAX2)
  150 CONTINUE
      DO 160 I = 1, N(5)
         MAX2 = MAX0(NLPB(I), MAX2)
  160 CONTINUE
      DO 170 I = 1, N(7)
         MAX3 = MAX0(NSPR(I), MAX3)
         MAX4 = MAX0(NHPR(I), MAX4)
  170 CONTINUE
      IF (REMESH) MAX1 = MAX1 * 20
      MAXNL = (MAX2 * (MAX3 + (MAX4 * MAX3))) + 1
      MAXNP = (MAX1 * MAXNL) + 1
      MAXPRM = 1 + MAX4
      MAX3 = MAX3 + 1

C  GET INITIAL SPACE IN ARRAY A FOR PERIMETER GENERATION

C  K(1) = X ARRAY OF THE PERIMETER
C  K(2) = Y ARRAY OF THE PERIMETER
C  K(3) = NID ARRAY OF THE PERIMETER
C  K(4) = LINE LIST
C  K(5) = NO OF NODES PER SIDE LIST
C  K(6) = WORK ARRAY FOR M1 GENERATION

      CALL MDRSRV ('X', K(1), MAXNP)
      CALL MDRSRV ('Y', K(2), MAXNP)
      CALL MDRSRV ('NID', K(3), MAXNP)
      CALL MDRSRV ('LISTL', K(4), MAXNL)
      CALL MDRSRV ('NNPS', K(5), MAX3)
      CALL MDRSRV ('ANGLE', K(6), MAXNP)
      CALL MDRSRV ('MARKED', K(26), MAXNL * 3)
      CALL MDSTAT (NERR, MUSED)
      IF (NERR .GT. 0) THEN
         CALL MDEROR (6)
         STOP ' '
      END IF

C  LOOP THROUGH THE GROUPS/REGIONS AND BAR SETS IN THE BODY LIST
C  CHECK CONNECTIVITY AND CALCULATE THE DIMENSIONS NEEDED FOR MESHING
C  NO PERIMETER INFORMATION IS SAVED THIS TIME THROUGH

      REAL = .FALSE.
      COUNT = .TRUE.
      DO 210 I = 1, N(9)
         CALL LTSORT (MR, LINKR, ABS(IRPB(I)), IPNTR1, ADDLNK)
         CALL LTSORT (MS, LINKB, ABS(IRPB(I)), IPNTR2, ADDLNK)

C  CHECK A REGION OR GROUP

         IF ((IRPB(I) .GT. 0) .AND. (IRPB(I) .LE. N(22)) .AND.
     &      (IPNTR1 .GT. 0)) THEN
            IF (IRGFLG(IPNTR1) .LE. -1) THEN
               WRITE (*, 10000) IRPB(I)
               L = IPNTR1

               CALL CHKRGN (IA, L, MP, ML, MS, MR, MSC, N(24), IPOINT,
     &            COOR, IPBOUN, ILINE, LTYPE, NINT, FACTOR, LCON,
     &            ILBOUN, ISBOUN, ISIDE, NLPS, IFLINE, ILLIST, IREGN,
     &            NSPR, IFSIDE, ISLIST, NPPF, IFPB, LISTPB, NLPF, IFLB,
     &            LISTLB, NSPF, IFSB, LISTSB, IFHOLE, NHPR, IHLIST,
     &            LINKP, LINKL, LINKS, LINKR, LINKSC, LINKPB, LINKLB,
     &            LINKSB, RSIZE, SCHEME, DEFSCH, NPREGN, NPSBC, NPNODE,
     &            MAXNP, MAXNL, MAX3, A(K(1)), A(K(2)), IA(K(3)),
     &            IA(K(4)), IA(K(5)), A(K(6)), A(K(26)), MXND, MXNPER,
     &            MXNL, MAXNBC, MAXSBC, AMESUR, XNOLD, YNOLD, NXKOLD,
     &            MMPOLD, LINKEG, LISTEG, BMESUR, MLINK, NPROLD, NPNOLD,
     &            NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX,
     &            IDIVIS, SIZMIN, EMAX, EMIN, NOROOM, ERRCHK, ERR)
            ELSE IF (IRGFLG(IPNTR1) .GE. 1) THEN
               WRITE (*, 10010) IRPB(I)
               J1 = IFSIDE(IPNTR1)
               J2 = J1 + NSPR(IPNTR1) - 1
               IREGN(IPNTR1) = -ABS(IREGN(IPNTR1))
               MXRXG = MAX(MXRXG, NSPR(IPNTR1))
               DO 180 J = J1, J2
                  CALL LTSORT (MR, LINKR, ABS(ISLIST(J)), IPNTR1,
     &               ADDLNK)
                  IF (IPNTR1 .GT. 0) THEN
                     WRITE (*, 10020) ABS(ISLIST(J))
                     L = IPNTR1

                     CALL CHKRGN (IA, L, MP, ML, MS, MR, MSC, N(24),
     &                  IPOINT, COOR, IPBOUN, ILINE, LTYPE, NINT,
     &                  FACTOR, LCON, ILBOUN, ISBOUN, ISIDE, NLPS,
     &                  IFLINE, ILLIST, IREGN, NSPR, IFSIDE, ISLIST,
     &                  NPPF, IFPB, LISTPB, NLPF, IFLB, LISTLB, NSPF,
     &                  IFSB, LISTSB, IFHOLE, NHPR, IHLIST, LINKP,
     &                  LINKL, LINKS, LINKR, LINKSC, LINKPB, LINKLB,
     &                  LINKSB, RSIZE, SCHEME, DEFSCH, NPREGN, NPSBC,
     &                  NPNODE, MAXNP, MAXNL, MAX3, A(K(1)),  A(K(2)),
     &                  IA(K(3)), IA(K(4)), IA(K(5)), A(K(6)),
     &                  A(K(26)), MXND, MXNPER, MXNL, MAXNBC, MAXSBC,
     &                  AMESUR, XNOLD, YNOLD, NXKOLD, MMPOLD, LINKEG,
     &                  LISTEG, BMESUR, MLINK, NPROLD, NPNOLD, NPEOLD,
     &                  NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX,
     &                  IDIVIS, SIZMIN, EMAX, EMIN, NOROOM, ERRCHK, ERR)
                  END IF
  180          CONTINUE
            END IF

C  WRITE AN ERROR FOR THIS REGION IN THE BODY LIST

         ELSE IF (IRPB(I) .GT. 0) THEN
            WRITE (*, 10030) IRPB(I)
            CALL LTSORT (MR, LINKR, ABS(IRPB(I)), IPNTR, ADDLNK)
            ADDLNK = .TRUE.
            IMINUS = -IPNTR
            CALL LTSORT (MR, LINKR, ABS(IRPB(I)), IMINUS, ADDLNK)
            ADDLNK = .FALSE.

C  CHECK A BAR SET

         ELSE IF ((IRPB(I) .LT. 0) .AND. (ABS(IRPB(I)) .LE. N(21))
     &      .AND. (IPNTR2 .GT. 0)) THEN
            L = IPNTR2
            WRITE (*, 10040) ABS(IRPB(I))
            REAL = .FALSE.
            COUNT = .TRUE.
            TEST = .FALSE.
            NPER = 1
            KNBC = 1
            KSBC = 1
            DO 190 J = JFLINE(IPNTR2), JFLINE(IPNTR2) + NLPB(IPNTR2) - 1
               CALL LTSORT (ML, LINKL, JLLIST(J), KK, ADDLNK)
               CALL LTSORT (MP, LINKP, LCON(1, KK), IP1, ADDLNK)
               CALL LTSORT (MP, LINKP, LCON(2, KK), IP2, ADDLNK)
               IF (LCON(3, KK) .GT. 0) THEN
                  CALL LTSORT (MP, LINKP, LCON(3, KK), IP3, ADDLNK)
               ELSE IF (LCON(3, KK) .LT. 0) THEN
                  CALL LTSORT (MP, LINKP, -LCON(3, KK), IP3, ADDLNK)
                  IP3 = -IP3
               ELSE
                  IP3 = 0
               END IF

               CALL PLINE (MP, ML, MAXNP, 1, 1, IPOINT, COOR, LINKP,
     &            ILINE(KK), LTYPE(KK), NINT(KK), FACTOR(KK), IP1, IP2,
     &            IP3, A(K(1)), A(K(2)), IA(K(3)), IPBOUN(IP1),
     &            IPBOUN(IP2), ILBOUN(KK), ISBOUN(KK), LINKPB, NPPF,
     &            IFPB, LISTPB, LINKLB, NLPF, IFLB, LISTLB, LINKSB,
     &            NSPF, IFSB, LISTSB, IDUMMY, KNBC, KSBC, ERR, TEST,
     &            REAL, COUNT, NOROOM, AMESUR, XNOLD, YNOLD, NXKOLD,
     &            MMPOLD, LINKEG, LISTEG, BMESUR, MLINK, NPROLD,
     &            NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN,
     &            REYMAX, IDIVIS, SIZMIN, EMAX, EMIN, GRAPH, DXMAX)
               IF (ERR) THEN
                  WRITE (*, 10050) IBARST(L)
                  ADDLNK = .FALSE.
                  CALL LTSORT (MS, LINKB, IRPB(I), IPNTR, ADDLNK)
                  ADDLNK = .TRUE.
                  IMINUS = -IPNTR
                  CALL LTSORT (MS, LINKB, IRPB(I), IMINUS, ADDLNK)
                  ADDLNK = .FALSE.
                  GO TO 200
               END IF
               NPER = NPER + NINT(KK) + 1
  190       CONTINUE
            IBARST(L) = -IBARST(L)

C  WHEN CHECKING THE MAXIMUMS - ADD ENOUGH FOR ONE MORE INTERVAL
C  ON THE LINE AS THIS LINE MAY BE INCREMENTED BY ONE IF THE
C  PERIMETER IS ODD

            MAXNBC = MAX0(MAXNBC, KNBC + 3)
            MAXSBC = MAX0(MAXSBC, KSBC + 3)
            MXND = MAX0(MXND, NPER)
            MXNPER = MAX0(MXNPER, NPER + 2)

C  WRITE AN ERROR FOR THIS BAR SET IN THE BODY LIST

         ELSE
            WRITE (*, 10060) ABS(IRPB(I))
            ADDLNK = .FALSE.
            CALL LTSORT (MS, LINKB, IRPB(I), IPNTR, ADDLNK)
            ADDLNK = .TRUE.
            IMINUS = -IPNTR
            CALL LTSORT (MS, LINKB, IRPB(I), IMINUS, ADDLNK)
            ADDLNK = .FALSE.
         END IF
  200    CONTINUE
  210 CONTINUE

C  RESET ALL USED POINTS AND LINES

      DO 220 I = 1, N(1)
         IPOINT(I) = ABS(IPOINT(I))
  220 CONTINUE
      DO 230 I = 1, N(2)
         NINT(I) = ABS(NINT(I))
  230 CONTINUE

C  RELEASE THE OLD ARRAYS, AND THEN
C  DIMENSION BASED ON THE MAXIMUMS CALCULATED

      CALL MDDEL ('X')
      CALL MDDEL ('Y')
      CALL MDDEL ('NID')
      CALL MDDEL ('LISTL')
      CALL MDDEL ('NNPS')
      CALL MDDEL ('ANGLE')
      CALL MDDEL ('MARKED')
      CALL MDSTAT (NERR, MUSED)
      IF (NERR .GT. 0) THEN
         CALL MDEROR (6)
         STOP ' '
      END IF

C  K(1) = X ARRAY OF THE PERIMETER
C  K(2) = Y ARRAY OF THE PERIMETER
C  K(3) = NID ARRAY(S) OF THE PERIMETER(S) [HOLES CAUSE MULTIPLE PERIMS]
C  K(4) = LINE LIST
C  K(5) = NO OF NODES PER SIDE LIST
C  K(6) = WORK ARRAY FOR M1 GENERATION
C  K(7) = XN (GLOBAL NODAL X VALUES)
C  K(8) = YN (GLOBAL NODAL Y VALUES)
C  K(9) = NUID (GLOBAL NODE UNIQUE ID'S)
C  K(10) = LXK (LINES PER ELEMENT)
C  K(11) = KXL (ELEMENTS PER LINE)
C  K(12) = NXL (NODES PER LINE)
C  K(13) = LXN (LINES PER NODE)
C  K(14) = LSTNBC (LIST OF NODAL BOUNDARY FLAGS AND NODES)
C  K(15) = LSTSBC (LIST OF SIDE BOUNDARY FLAGS AND NODES)
C  K(16) = XSUB ARRAY OF THE PERIMETER OF A SUBREGION
C  K(17) = YSUB ARRAY OF THE PERIMETER OF A SUBREGION
C  K(18) = NIDSUB ARRAY OF THE PERIMETER OF A SUBREGION
C  K(19) = NODE NUMBERS AROUND HOLE VOID
C  K(20) = NUMBER OF NODES ON PERIMETERS (REGION + HOLE)
C  K(21) = INDEX ARRAY FOR COMBINING SUB-REGIONS AND REGIONS
C  K(22) = FANGLE ARRAY FOR INTERIOR ANGLE IN FILL ROUTINES
C  K(23) = BNSIZE ARRAY FOR SIZE DIFFERENTIAL IN FILL ROUTINES
C  K(24) = LNODES ARRAY FOR CONNECTIVITY OF THE INSIDE PERIMETER
C          NODES IN FILL ROUTINES
C  NOTE: LINES IN THIS CONTEXT REFERS TO CONNECTIONS OF ELEMENT NODES

C  MAKE ROOM IN LINE LIST FOR HOLES

      MXND = INT(MXND * MXRXG * 1.2)
      MXNL = MXNL + ( (MXRXG + MAX4) * MAX2 * MAX3)
      MLN = 8
      MAXNB = MXNPER * MAXPRM
      CALL MDRSRV ('X', K(1), MXNPER)
      CALL MDRSRV ('Y', K(2), MXNPER)
      CALL MDRSRV ('NID', K(3), MXNPER * MAXPRM)
      CALL MDRSRV ('LISTL', K(4), MXNL)
      CALL MDRSRV ('NNPS', K(5), MAX3)
      CALL MDRSRV ('ANGLE', K(6), MXNPER)
      CALL MDRSRV ('XN', K(7), MXND)
      CALL MDRSRV ('YN', K(8), MXND)
      CALL MDRSRV ('NUID', K(9), MXND)
      CALL MDRSRV ('LXK', K(10), MXND*4)
      CALL MDRSRV ('KXL', K(11), MXND*6)
      CALL MDRSRV ('NXL', K(12), MXND*6)
      CALL MDRSRV ('LXN', K(13), MXND*4)
      CALL MDRSRV ('LSTNBC', K(14), MAXNBC)
      CALL MDRSRV ('LSTSBC', K(15), 2*MAXSBC)
      CALL MDRSRV ('XSUB', K(16), MXNPER)
      CALL MDRSRV ('YSUB', K(17), MXNPER)
      CALL MDRSRV ('NIDSUB', K(18), MXNPER)
      CALL MDRSRV ('NXH', K(19), MXND)
      CALL MDRSRV ('NPERIM', K(20), MAXPRM)
      CALL MDRSRV ('INDX', K(21), MXND)
      CALL MDRSRV ('FANGLE', K(22), MXND)
      CALL MDRSRV ('BNSIZE', K(23), MXND * 2)
      CALL MDRSRV ('LNODES', K(24), MXND * MLN)
      CALL MDRSRV ('PRLINK', K(25), MAXPRM * 3)
      CALL MDRSRV ('MARKED', K(26), MXNL * 3)
      CALL MDRSRV ('IPTPER', K(27), MAXPRM)
      CALL MDRSRV ('NUMPER', K(28), MAXPRM)
      CALL MDRSRV ('LPERIM', K(29), MAXNB)
      CALL MDRSRV ('ZN', K(30), MXND)
      CALL MDSTAT (NERR, MUSED)
      IF (NERR .GT. 0) THEN
         CALL MDEROR (6)
         STOP ' '
      END IF

C  SET UP THE LOOP FOR PROCESSING GROUPS

      IF (LGROUP) THEN
  240    CONTINUE
         IF (STEP .AND. (N(22) .GT. 0)) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('STEP PROCESS GROUPS I1 THROUGH I2')
            IF (ICOM .GT. JCOM) THEN
               CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN,
     &            IIN, RIN)
               ICOM = 1
            END IF
            CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2,
     &         IFOUND)
            IF (IFOUND .GT. 0) THEN
               CALL CHECK (I1, I2, N(22))
            ELSE
               GO TO 320
            END IF
         ELSE
            I1 = 1
            I2 = N(22)
         END IF

C  BEGIN PROCESSING GROUPS

         REAL = .TRUE.
         COUNT = .FALSE.
         DO 310 IGRP = I1, I2
  250       CONTINUE
            CALL LTSORT (MR, LINKR, IGRP, IGPNTR, ADDLNK)
            IF ((IGPNTR .GT. 0) .AND. (IRGFLG(IGPNTR) .GE. 1)
     &         .AND. (IREGN(IGPNTR) .LT. 0)) THEN
               WRITE (*, 10070) IGRP
               J1 = IFSIDE(IGPNTR)
               J2 = J1 + NSPR(IGPNTR) - 1
               NNN = 0
               KKK = 0
               LLL = 0
               NNNOLD = 0
               KKKOLD = 0
               LLLOLD = 0
               DO 270 J = J1, J2
                  CALL LTSORT (MR, LINKR, ABS(ISLIST(J)), IPNTR2,
     &               ADDLNK)
                  IF ((IPNTR2 .GT. 0) .AND. (IREGN(IPNTR2) .LT. 0)) THEN
                     L = IPNTR2
                     NOROOM = .FALSE.
                     CALL MESAGE (' ')
                     WRITE (*, 10080) ABS(IREGN(L))

C  CALCULATE THE PERIMETER OF THE REGION

  260                CONTINUE
                     NPRM = 1
                     JJHOLE = 0
                     KNBC = 0

                     CALL PERIM (MP, ML, MS, NSPR(L), MXNL, MXNPER,
     &                  MAXNBC, MAXSBC, KNBC, KSBC, ABS (IREGN(L)),
     &                  IPOINT, COOR, IPBOUN, ILINE, LTYPE, NINT,
     &                  FACTOR, LCON, ILBOUN, ISBOUN, ISIDE, NLPS,
     &                  IFLINE, ILLIST, ISLIST(IFSIDE(L)), NPPF, IFPB,
     &                  LISTPB, NLPF, IFLB, LISTLB, NSPF, IFSB, LISTSB,
     &                  LINKP, LINKL, LINKS, LINKPB, LINKLB, LINKSB,
     &                  A(K(1)), A(K(2)), IA(K(3)), NPER, IA(K(4)), NL,
     &                  IA(K(14)), IA(K(26)), EVEN, REAL, ERR, CCW,
     &                  COUNT, NOROOM, AMESUR, XNOLD, YNOLD, NXKOLD,
     &                  MMPOLD, LINKEG, LISTEG, BMESUR, MLINK, NPROLD,
     &                  NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX,
     &                  REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN)

C  GET THE REGION SCHEME

                     CALL LTSORT (MR, LINKSC, ABS(IREGN(L)), IPNTR,
     &                  ADDLNK)
                     CALL RGNSCH (MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN,
     &                  STEP, IREGN(L), IPNTR, N(24), MSC, SCHEME,
     &                  DEFSCH, SCHSTR, LENSCH, NPER, PENTAG, TRIANG,
     &                  TRNSIT, HALFC, FILL, ICODE, REMESH)
                     IF (ICODE .EQ. IEXIT) THEN
                        GO TO 270
                     ELSE IF (ICODE .EQ. IOVER) THEN
                        GO TO 260
                     ELSE IF (ICODE .EQ. IQUIT) THEN
                        GO TO 270

C  GENERATE INITIAL GRID

C  CALCULATE A "TRANSITION" MAPPED MESH

                     ELSE IF (TRNSIT) THEN
                        CALL BMSCHM (NPER, KKK, LLL, NNN, ML, MS,
     &                     NSPR(L), ISLIST(IFSIDE(L)), NINT, IFLINE,
     &                     NLPS, ILLIST, LINKL, LINKS, MXNPER, MAXPRM,
     &                     MAX3, MXND, A(K(1)), A(K(2)), IA(K(3)),
     &                     IA(K(5)), A(K(6)), A(K(7)), A(K(8)),
     &                     IA(K(9)), IA(K(10)), IA(K(11)), IA(K(12)),
     &                     IA(K(13)), A(K(16)), A(K(17)), IA(K(18)),
     &                     IA(K(21)), IAVAIL, NAVAIL, CCW, HALFC, ERR)

C  CALCULATE A "TRIANGULAR" MAPPED MESH

                     ELSE IF (TRIANG) THEN
                        CALL TMSCHM (NPER, KKK, LLL, NNN, ML, MS,
     &                     NSPR(L), ISLIST(IFSIDE(L)), NINT, IFLINE,
     &                     NLPS, ILLIST, LINKL, LINKS, MXNPER, MAXPRM,
     &                     MAX3, MXND, A(K(1)), A(K(2)), IA(K(3)),
     &                     IA(K(5)), A(K(6)), A(K(7)), A(K(8)),
     &                     IA(K(9)), IA(K(10)), IA(K(11)), IA(K(12)),
     &                     IA(K(13)), A(K(16)), A(K(17)), IA(K(18)),
     &                     IA(K(21)), IAVAIL, NAVAIL, CCW, ERR)

C  CALCULATE A "PENTAGON" MAPPED MESH

                     ELSE IF (PENTAG) THEN
                        CALL UMSCHM (IA, NPER, KKK, LLL, NNN, ML, MS,
     &                     NSPR(L), ISLIST(IFSIDE(L)), NINT, IFLINE,
     &                     NLPS, ILLIST, LINKL, LINKS, MXNPER, MAXPRM,
     &                     MAX3, MXND, A(K(1)), A(K(2)), IA(K(3)),
     &                     IA(K(5)), A(K(6)), A(K(7)), A(K(8)),
     &                     IA(K(9)), IA(K(10)), IA(K(11)), IA(K(12)),
     &                     IA(K(13)), A(K(16)), A(K(17)), IA(K(18)),
     &                     IA(K(21)), IAVAIL, NAVAIL, CCW, ERR)

C  USE THE PAVING TECHNIQUE TO FILL THE INITIAL REGION

                     ELSE IF (FILL) THEN

                        CALL PMSCHM (NPER, NPRM, MXND, MLN, MP, ML, MS,
     &                     MR, NL, MXNL, MXNPER, MAXPRM, MAXNB, MAXNBC,
     &                     MAXSBC, KNBC, KSBC, KNUM, IPOINT, COOR,
     &                     IPBOUN, ILINE, LTYPE, NINT, FACTOR, LCON,
     &                     ILBOUN, ISBOUN, ISIDE, NLPS, IFLINE, ILLIST,
     &                     ISLIST, IREGN, NPPF, IFPB, LISTPB, NLPF,
     &                     IFLB, LISTLB, NSPF, IFSB, LISTSB, LINKP,
     &                     LINKL, LINKS, LINKR, LINKPB, LINKLB, LINKSB,
     &                     NSPR, IFSIDE, RSIZE, IFHOLE, NHPR, IHLIST,
     &                     A(K(1)), A(K(2)), IA(K(3)), IA(K(4)),
     &                     A(K(7)), A(K(8)), A(K(30)), IA(K(9)),
     &                     IA(K(10)), IA(K(11)), IA(K(12)), IA(K(13)),
     &                     IA(K(14)), IA(K(20)), A(K(22)), A(K(23)),
     &                     IA(K(24)), IA(K(25)), IA(K(26)), IA(K(27)),
     &                     IA(K(28)), IA(K(29)), KKK, NNN, LLL, IAVAIL,
     &                     NAVAIL, DEV1, ABS(IREGN(L)), L, BATCH,NOROOM,
     &                     ERR, AMESUR, XNOLD, YNOLD, NXKOLD, MMPOLD,
     &                     LINKEG, LISTEG, BMESUR, MLINK, NPROLD,
     &                     NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX,
     &                     REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN)
                        IF (NOROOM) THEN
                           MXND = INT(MXND*1.5 + 1)
                           MAXNBC = MAX0 (MAXNBC, KNBC)
                           MAXSBC = MAX0 (MAXSBC, KSBC)
                           KSBC = KKSBC
                           CALL MDLONG ('XN', K(7), MXND)
                           CALL MDLONG ('YN', K(8), MXND)
                           CALL MDLONG ('ZN', K(30), MXND)
                           CALL MDLONG ('NUID', K(9), MXND)
                           CALL MDLONG ('LXK', K(10), MXND*4)
                           CALL MDLONG ('KXL', K(11), MXND*6)
                           CALL MDLONG ('NXL', K(12), MXND*6)
                           CALL MDLONG ('LXN', K(13), MXND*4)
                           CALL MDLONG ('LSTNBC', K(14), MAXNBC)
                           CALL MDLONG ('LSTSBC', K(15), MAXSBC)
                           CALL MDLONG ('NXH', K(19), MXND)
                           CALL MDLONG ('INDX', K(21), MXND)
                           CALL MDLONG ('FANGLE', K(22), MXND)
                           CALL MDLONG ('BNSIZE', K(23), MXND * 2)
                           CALL MDLONG ('LNODES', K(24), MXND * MLN)
                           CALL MDSTAT (NERR, MUSED)
                           IF (NERR .GT. 0) THEN
                              CALL MDEROR (6)
                              STOP ' '
                           END IF
                           CALL MESAGE
     &                        ('REDIMENSIONING NEEDED - PLEASE WAIT')
                           IF (STEP) THEN
                              CALL MESAGE
     &                           ('CURRENT PROCESSING SCHEME IS SAVED')
                           ELSE
                              CALL MESAGE
     &                           ('CURRENT SCHEME WILL BE REPEATED')
                           END IF
                           GO TO 260
                        END IF

C  PROCESS A "NORMAL" REGION

                     ELSE
                        CALL MMSCHM (NPER, KKK, LLL, NNN, ML, MS,
     &                     NSPR(L), ISLIST(IFSIDE(L)), NINT, IFLINE,
     &                     NLPS, ILLIST, LINKL, LINKS, MXNPER, MAXPRM,
     &                     MAX3, MXND, A(K(1)), A(K(2)), IA(K(3)),
     &                     IA(K(5)), A(K(6)), A(K(7)), A(K(8)),
     &                     IA(K(9)), IA(K(10)), IA(K(11)), IA(K(12)),
     &                     IA(K(13)), IAVAIL, NAVAIL, CCW, REAL, SCHSTR,
     &                     M1, ERR)
                     END IF

C  FLAG THE REGION IF AN ERROR HAS OCCURRED

                     IF (ERR) THEN
                        CALL MESAGE ('ERROR IN INITIAL GRID GENERATION')
                        CALL MESAGE ('** REGION PROCESSING ABORTED **')
                        CALL PLTBEL
                        CALL PLTFLU
                        CALL MESAGE (' ')
                        IREGN(L) = ABS(IREGN(L))
                        IF (ISLIST(J) .EQ. IREGN(L)) THEN
                           ADDLNK = .TRUE.
                           IMINUS = -L
                           CALL LTSORT (MR, LINKR, IREGN(L), IMINUS,
     &                        ADDLNK)
                           ADDLNK = .FALSE.
                        END IF
                        GO TO 270
                     END IF

C  BEGIN FULL SCHEME CONTROL FOR A GROUP SUB-REGION

                     RECT = .NOT.(PENTAG .OR. TRIANG .OR.
     &                  TRNSIT .OR. FILL)
                     IF (STEP) CALL MINMAX_FQ (MXNPER, NPER, A(K(1)),
     &                  A(K(2)), XMIN, XMAX, YMIN, YMAX)

                     CALL PSCHEM (MP, ML, MS, MR, N, IPOINT, COOR,
     &                  IPBOUN, ILINE, LTYPE, NINT, FACTOR, LCON,
     &                  ILBOUN, ISBOUN, ISIDE, NLPS, IFLINE, ILLIST,
     &                  IREGN, NSPR, IFSIDE, ISLIST, NPPF, IFPB, LISTPB,
     &                  NLPF, IFLB, LISTLB, NSPF, IFSB, LISTSB, LINKP,
     &                  LINKL, LINKS, LINKR, LINKSC, LINKPB, LINKLB,
     &                  LINKSB, IFHOLE, NHPR, IHLIST, MAXNBC, KNBC,
     &                  MAXSBC, KSBC, MXND, NNN, NNNOLD, KKK, KKKOLD,
     &                  LLL, A(K(1)), A(K(2)), IA(K(3)), IA(K(4)),
     &                  A(K(7)), A(K(8)), IA(K(9)), IA(K(10)),
     &                  IA(K(11)), IA(K(12)), IA(K(13)), IA(K(14)),
     &                  IA(K(19)), IA(K(20)), IA(K(26)), IAVAIL, NAVAIL,
     &                  MXNL, MXNPER, NPER, MAXPRM, NPRM, MSC, ISCHM,
     &                  SCHEME, SCHSTR, RECT, M1, INSIDE, JJHOLE, KKSBC,
     &                  DEV1, EIGHT, NINE, STEP, L, NL, MCOM, CIN, IIN,
     &                  RIN, KIN, ICOM, JCOM, XMIN, XMAX, YMIN, YMAX,
     &                  ICODE, NOROOM, ERR, AMESUR, XNOLD, YNOLD,NXKOLD,
     &                  MMPOLD, LINKEG, LISTEG, BMESUR, MLINK, NPROLD,
     &                  NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX,
     &                  REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN)
                     IF (NOROOM) THEN
                        MXND = INT(MXND*1.2 + 1)
                        MAXNBC = MAX0 (MAXNBC, KNBC)
                        MAXSBC = MAX0 (MAXSBC, KSBC)
                        KSBC = KKSBC
                        CALL MDLONG ('XN', K(7), MXND)
                        CALL MDLONG ('YN', K(8), MXND)
                        CALL MDLONG ('ZN', K(30), MXND)
                        CALL MDLONG ('NUID', K(9), MXND)
                        CALL MDLONG ('LXK', K(10), MXND*4)
                        CALL MDLONG ('KXL', K(11), MXND*6)
                        CALL MDLONG ('NXL', K(12), MXND*6)
                        CALL MDLONG ('LXN', K(13), MXND*4)
                        CALL MDLONG ('LSTNBC', K(14), MAXNBC)
                        CALL MDLONG ('LSTSBC', K(15), MAXSBC)
                        CALL MDLONG ('NXH', K(19), MXND)
                        CALL MDLONG ('INDX', K(21), MXND)
                        CALL MDLONG ('FANGLE', K(22), MXND)
                        CALL MDLONG ('BNSIZE', K(23), MXND * 2)
                        CALL MDLONG ('LNODES', K(24), MXND * MLN)
                        CALL MDSTAT (NERR, MUSED)
                        IF (NERR .GT. 0) THEN
                           CALL MDEROR (6)
                           STOP ' '
                        END IF
                        CALL MESAGE
     &                     ('REDIMENSIONING NEEDED - PLEASE WAIT')
                        IF (STEP) THEN
                           CALL MESAGE
     &                        ('CURRENT PROCESSING SCHEME IS SAVED')
                        ELSE
                           CALL MESAGE
     &                        ('CURRENT SCHEME WILL BE REPEATED')
                        END IF
                        GO TO 260
                     ELSE IF (ERR) THEN
                        IF (STEP) THEN
                           CALL INTRUP
     &                        ('WOULD YOU LIKE TO REPROCESS REGION',
     &                        IANS, MCOM, ICOM, JCOM, CIN, IIN, RIN,
     &                        KIN)
                           NNN = NNNOLD
                           KKK = KKKOLD
                           LLL = LLLOLD
                           IF (IANS) GO TO 260
                        END IF
                        CALL MESAGE ('REGION PROCESSING ABORTED')
                        GO TO 270
                     ELSE IF (ICODE .EQ. IEXIT) THEN
                        FINAL = J .EQ. J2
                        IF (J .GT. J1) THEN
                           CALL FIXSUB (MXND, NNNOLD, NNN, LLLOLD, LLL,
     &                        KKKOLD, KKK, A(K(7)), A(K(8)), IA(K(9)),
     &                        IA(K(10)), IA(K(11)), IA(K(12)),
     &                        IA(K(13)), IA(K(21)), IAVAIL, NAVAIL,
     &                        FINAL)
                        END IF
                        NNNOLD = NNN
                        LLLOLD = LLL
                        KKKOLD = KKK
                        IF (FINAL) CALL FXNUID (NSPR(IGPNTR),
     &                     ISLIST(IFSIDE(IGPNTR)), MR, MS, ML, NSPR,
     &                     ILINE, ISIDE, NLPS, IFLINE, ILLIST, LCON,
     &                     ISLIST, IFSIDE, LINKR, LINKS, LINKL, NNN,
     &                     MXNL, MXND, IA(K(4)), IA(K(9)), IA(K(12)),
     &                     IA(K(13)), IA(K(21)), NOROOM, ERR)
                        IF (ERR) THEN
                           CALL MESAGE
     &                        ('GROUP SCHEME PROCESSING NOT POSSIBLE')
                           CALL MESAGE ('GROUP PROCESSING ABORTED')
                           GO TO 300
                        END IF
                     ELSE IF (ICODE .EQ. IOVER) THEN
                        NNN = NNNOLD
                        KKK = KKKOLD
                        LLL = LLLOLD
                        GO TO 260
                     ELSE IF (ICODE .EQ. IQUIT) THEN
                        NNN = NNNOLD
                        KKK = KKKOLD
                        LLL = LLLOLD
                        GO TO 270
                     END IF
                  END IF
  270          CONTINUE

C  BEGIN FULL SCHEME CONTROL FOR A GROUP REGION

               NNNOLD = 0
               KKKOLD = 0
               RECT = .FALSE.
               CALL MESAGE (' ')
               CALL MESAGE ('GROUP SCHEME PROCESSING BEGUN')
               CALL LTSORT (MR, LINKSC, ABS(IREGN(IGPNTR)), IPNTR,
     &            ADDLNK)
               CALL RGNSCH (MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN, STEP,
     &            IREGN(IGPNTR), IPNTR, N(24), MSC, SCHEME, DEFSCH,
     &            SCHSTR, LENSCH, NPER, PENTAG, TRIANG, TRNSIT, HALFC,
     &            FILL, ICODE, REMESH)

               IF (ICODE .EQ. IEXIT) THEN
                  CALL CHKKXL (MXND, IA(K(10)), IA(K(11)), LLL, ERR)
                  IF (ERR) THEN
                     CALL MESAGE ('ERROR IN CHECK OF KXL ARRAY')
                     IF (STEP) THEN
                        CALL INTRUP
     &                     ('WOULD YOU LIKE TO REPROCESS GROUP',
     &                     IANS, MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN)
                        IF (IANS) GO TO 250
                     END IF
                     CALL MESAGE ('GROUP PROCESSING ABORTED')
                     GO TO 300
                  END IF
                  BAR = .FALSE.
                  KSBC = 0
                  CALL GETSBC (MXND, MXNPER, NPER, NL, ML, MAXSBC,
     &               MAXPRM, NPRM, IA(K(3)), IA(K(4)), A(K(7)), A(K(8)),
     &               IA(K(9)), IA(K(10)), IA(K(11)), IA(K(12)),
     &               IA(K(15)), IA(K(20)), KSBC, LCON, ISBOUN, LINKL,
     &               NSPF, IFSB, LISTSB, LINKSB, LLL, BAR, ERR)
                  IF (ERR) THEN
                     CALL MESAGE ('ERROR IN SORTING SIDE BOUNDARIES')
                     IF (STEP) THEN
                        CALL INTRUP
     &                     ('WOULD YOU LIKE TO REPROCESS GROUP',
     &                     IANS, MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN)
                        IF (IANS) GO TO 250
                     END IF
                     CALL MESAGE ('GROUP PROCESSING ABORTED')
                     GO TO 300
                  END IF
                  CALL MKUSED (MXNL, MP, ML, IA(K(4)), IPOINT, NINT,
     &               LINKP, LINKL, LCON, NL)
                  CALL SAVREG (MXND, MAXNBC, MAXSBC, A(K(7)), A(K(8)),
     &               IA(K(9)), IA(K(10)), IA(K(12)), IA(K(13)),
     &               IA(K(14)), IA(K(15)), KNBC, KSBC, NNN, KKK, IGRP,
     &               IUNIT, BAR, M1)
                  DO 280 J = J1, J2
                     CALL LTSORT (MR, LINKR, ABS(ISLIST(J)), IPNTR2,
     &                  ADDLNK)
                     IREGN(IPNTR2) = ABS(IREGN(IPNTR2))
  280             CONTINUE
                  NPREGN = NPREGN + 1
                  NPELEM = NPELEM + KKK
                  NPNODE = NPNODE + NNN
                  IF (EIGHT .OR. NINE) NPNODE = NPNODE + LLL
                  IF (NINE) NPNODE = NPNODE + KKK
                  NPNBC = NPNBC + KNBC
                  IF (THREE .OR. EIGHT .OR. NINE)
     &               NPNBC = NPNBC + KNBC
                  NPSBC = NPSBC + KSBC
                  IF (THREE .OR. EIGHT .OR. NINE)
     &               NPSBC = NPSBC + KSBC
                  MAXKXN = MAXKXN + LLL
                  IREGN(IGPNTR) = ABS(IREGN(IGPNTR))
                  WRITE (*, 10100) IREGN(IGPNTR)
                  GO TO 300
               ELSE IF (ICODE .EQ. IOVER) THEN
                  GO TO 250
               ELSE IF (ICODE .EQ. IQUIT) THEN
                  GO TO 300
               END IF

               IF (STEP) CALL MINMAX_FQ (MXND, NNN, A(K(7)), A(K(8)),
     &            XMIN, XMAX, YMIN, YMAX)

               CALL PSCHEM (MP, ML, MS, MR, N, IPOINT, COOR, IPBOUN,
     &            ILINE, LTYPE, NINT, FACTOR, LCON, ILBOUN, ISBOUN,
     &            ISIDE, NLPS, IFLINE, ILLIST, IREGN, NSPR, IFSIDE,
     &            ISLIST, NPPF, IFPB, LISTPB, NLPF, IFLB, LISTLB, NSPF,
     &            IFSB, LISTSB, LINKP, LINKL, LINKS, LINKR, LINKSC,
     &            LINKPB, LINKLB, LINKSB, IFHOLE, NHPR, IHLIST, MAXNBC,
     &            KNBC, MAXSBC, KSBC, MXND, NNN, NNNOLD, KKK, KKKOLD,
     &            LLL, A(K(1)), A(K(2)), IA(K(3)), IA(K(4)), A(K(7)),
     &            A(K(8)), IA(K(9)), IA(K(10)), IA(K(11)), IA(K(12)),
     &            IA(K(13)), IA(K(14)), IA(K(19)), IA(K(20)), IA(K(26)),
     &            IAVAIL, NAVAIL, MXNL, MXNPER, NPER, MAXPRM, NPRM, MSC,
     &            ISCHM, SCHEME, SCHSTR, RECT, M1, INSIDE, JJHOLE,
     &            KKSBC, DEV1, EIGHT, NINE, STEP, IGPNTR, NL, MCOM, CIN,
     &            IIN, RIN, KIN, ICOM, JCOM, XMIN, XMAX, YMIN, YMAX,
     &            ICODE, NOROOM, ERR, AMESUR, XNOLD, YNOLD, NXKOLD,
     &            MMPOLD, LINKEG, LISTEG, BMESUR, MLINK, NPROLD, NPNOLD,
     &            NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX,
     &            IDIVIS, SIZMIN, EMAX, EMIN)
               IF (NOROOM) THEN
                  MXND = INT(MXND*1.2 + 1)
                  MAXNBC = MAX0 (MAXNBC, KNBC)
                  MAXSBC = MAX0 (MAXSBC, KSBC)
                  KSBC = KKSBC
                  CALL MDLONG ('XN', K(7), MXND)
                  CALL MDLONG ('YN', K(8), MXND)
                  CALL MDLONG ('ZN', K(30), MXND)
                  CALL MDLONG ('NUID', K(9), MXND)
                  CALL MDLONG ('LXK', K(10), MXND*4)
                  CALL MDLONG ('KXL', K(11), MXND*6)
                  CALL MDLONG ('NXL', K(12), MXND*6)
                  CALL MDLONG ('LXN', K(13), MXND*4)
                  CALL MDLONG ('LSTNBC', K(14), MAXNBC)
                  CALL MDLONG ('LSTSBC', K(15), MAXSBC)
                  CALL MDLONG ('NXH', K(19), MXND)
                  CALL MDLONG ('INDX', K(21), MXND)
                  CALL MDLONG ('FANGLE', K(22), MXND)
                  CALL MDLONG ('BNSIZE', K(23), MXND * 2)
                  CALL MDLONG ('LNODES', K(24), MXND * MLN)
                  CALL MDSTAT (NERR, MUSED)
                  IF (NERR .GT. 0) THEN
                     CALL MDEROR (6)
                     STOP ' '
                  END IF
                  CALL MESAGE
     &               ('REDIMENSIONING NEEDED - PLEASE WAIT')
                  IF (STEP) THEN
                     CALL MESAGE ('CURRENT PROCESSING SCHEME IS SAVED')
                  ELSE
                     CALL MESAGE ('CURRENT SCHEME WILL BE REPEATED')
                  END IF
                  GO TO 250
               ELSE IF (ERR) THEN
                  IF (STEP) THEN
                     CALL INTRUP ('WOULD YOU LIKE TO REPROCESS GROUP',
     &                  IANS, MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN)
                     IF (IANS) GO TO 250
                  END IF
                  CALL MESAGE ('GROUP PROCESSING ABORTED')
                  GO TO 300
               ELSE IF (ICODE .EQ. IEXIT) THEN
                  CALL CHKKXL (MXND, IA(K(10)), IA(K(11)), LLL, ERR)
                  IF (ERR) THEN
                     CALL MESAGE ('ERROR IN CHECK OF KXL ARRAY')
                     IF (STEP) THEN
                        CALL INTRUP
     &                     ('WOULD YOU LIKE TO REPROCESS GROUP',
     &                     IANS, MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN)
                        IF (IANS) GO TO 250
                     END IF
                     CALL MESAGE ('GROUP PROCESSING ABORTED')
                     GO TO 300
                  END IF
                  BAR = .FALSE.
                  KSBC = 0
                  CALL GETSBC (MXND, MXNPER, NPER, NL, ML, MAXSBC,
     &               MAXPRM, NPRM, IA(K(3)), IA(K(4)), A(K(7)), A(K(8)),
     &               IA(K(9)), IA(K(10)), IA(K(11)), IA(K(12)),
     &               IA(K(15)), IA(K(20)), KSBC, LCON, ISBOUN, LINKL,
     &               NSPF, IFSB, LISTSB, LINKSB, LLL, BAR, ERR)
                  IF (ERR) THEN
                     CALL MESAGE ('ERROR IN SORTING SIDE BOUNDARIES')
                     IF (STEP) THEN
                        CALL INTRUP
     &                     ('WOULD YOU LIKE TO REPROCESS GROUP',
     &                     IANS, MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN)
                        IF (IANS) GO TO 250
                     END IF
                     CALL MESAGE ('GROUP PROCESSING ABORTED')
                     GO TO 300
                  END IF
                  CALL MKUSED (MXNL, MP, ML, IA(K(4)), IPOINT, NINT,
     &               LINKP, LINKL, LCON, NL)
                  CALL SAVREG (MXND, MAXNBC, MAXSBC, A(K(7)), A(K(8)),
     &               IA(K(9)), IA(K(10)), IA(K(12)), IA(K(13)),
     &               IA(K(14)), IA(K(15)), KNBC, KSBC, NNN, KKK, IGRP,
     &               IUNIT, BAR, M1)
                  DO 290 J = J1, J2
                     CALL LTSORT (MR, LINKR, ABS(ISLIST(J)), IPNTR2,
     &                  ADDLNK)
                     IREGN(IPNTR2) = ABS(IREGN(IPNTR2))
  290             CONTINUE
                  NPREGN = NPREGN + 1
                  NPELEM = NPELEM + KKK
                  NPNODE = NPNODE + NNN
                  IF (EIGHT .OR. NINE) NPNODE = NPNODE + LLL
                  IF (NINE) NPNODE = NPNODE + KKK
                  NPNBC = NPNBC + KNBC
                  IF (THREE .OR. EIGHT .OR. NINE)
     &               NPNBC = NPNBC + KNBC
                  NPSBC = NPSBC + KSBC
                  IF (THREE .OR. EIGHT .OR. NINE)
     &               NPSBC = NPSBC + KSBC
                  MAXKXN = MAXKXN + LLL
                  IREGN(IGPNTR) = ABS(IREGN(IGPNTR))
                  WRITE (*, 10100) IREGN(IGPNTR)
               ELSE IF (ICODE .EQ. IOVER) THEN
                  GO TO 250
               ELSE IF (ICODE .EQ. IQUIT) THEN
                  GO TO 300
               END IF
            END IF
  300       CONTINUE
  310    CONTINUE

C  END OF THIS SET OF GROUPS
C     IF STEPPING THROUGH, SEE IF ANY MORE GROUPS
C        ARE TO BE PROCESSED

         IF (STEP) THEN
            CALL INTRUP ('PROCESS ADDITIONAL GROUPS', IANS, MCOM,
     &         ICOM, JCOM, CIN, IIN, RIN, KIN)
            IF (IANS) GO TO 240
         END IF
      END IF

C  SET UP THE LOOP FOR PROCESSING REGIONS

  320 CONTINUE
      IF (STEP .AND. (N(22) .GT. 0)) THEN
         CALL MESAGE (' ')
         CALL MESAGE ('STEP PROCESS REGIONS I1 THROUGH I2')
         IF (ICOM .GT. JCOM) THEN
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN, IIN,
     &         RIN)
            ICOM = 1
         END IF
         CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2, IFOUND)
         IF (IFOUND .GT. 0) THEN
            CALL CHECK (I1, I2, N(22))
         ELSE
            GO TO 370
         END IF
      ELSE
         I1 = 1
         I2 = N(22)
      END IF

C  BEGIN PROCESSING REGIONS

      REAL = .TRUE.
      COUNT = .FALSE.
      DO 350 I = I1, I2
         CALL LTSORT (MR, LINKR, I, L, ADDLNK)
         IF (L .GT. 0) THEN
         IF (IREGN(L) .LT. 0) THEN
            NOROOM = .FALSE.
            CALL MESAGE (' ')
            WRITE (*, 10090) ABS(IREGN(L))

C  CALCULATE THE PERIMETER OF THE REGION

  330       CONTINUE
            NNN = 0
            KKK = 0
            LLL = 0
            NPRM = 1
            JJHOLE = 0
            KNBC = 0
            CALL PERIM (MP, ML, MS, NSPR(L), MXNL, MXNPER, MAXNBC,
     &         MAXSBC, KNBC, KSBC, ABS (IREGN(L)), IPOINT, COOR, IPBOUN,
     &         ILINE, LTYPE, NINT, FACTOR, LCON, ILBOUN, ISBOUN, ISIDE,
     &         NLPS, IFLINE, ILLIST, ISLIST(IFSIDE(L)), NPPF, IFPB,
     &         LISTPB, NLPF, IFLB, LISTLB, NSPF, IFSB, LISTSB, LINKP,
     &         LINKL, LINKS, LINKPB, LINKLB, LINKSB, A(K(1)), A(K(2)),
     &         IA(K(3)), NPER, IA(K(4)), NL, IA(K(14)), IA(K(26)), EVEN,
     &         REAL, ERR, CCW, COUNT, NOROOM, AMESUR, XNOLD, YNOLD,
     &         NXKOLD, MMPOLD, LINKEG, LISTEG, BMESUR, MLINK, NPROLD,
     &         NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN,
     &         REYMAX, IDIVIS, SIZMIN, EMAX, EMIN)

C  GET THE REGION SCHEME

            CALL LTSORT (MR, LINKSC, ABS(IREGN(L)), IPNTR, ADDLNK)
            CALL RGNSCH (MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN, STEP,
     &         IREGN(L), IPNTR, N(24), MSC, SCHEME, DEFSCH, SCHSTR,
     &         LENSCH, NPER, PENTAG, TRIANG, TRNSIT, HALFC, FILL, ICODE,
     &         REMESH)
            IF (ICODE .EQ. IEXIT) THEN
               GO TO 360
            ELSE IF (ICODE .EQ. IOVER) THEN
               GO TO 330
            ELSE IF (ICODE .EQ. IQUIT) THEN
               GO TO 350

C  GENERATE INITIAL GRID

C  CALCULATE A "TRANSITION" MAPPED MESH

            ELSE IF (TRNSIT) THEN
               CALL BMSCHM (NPER, KKK, LLL, NNN, ML, MS, NSPR(L),
     &            ISLIST(IFSIDE(L)), NINT, IFLINE, NLPS, ILLIST, LINKL,
     &            LINKS, MXNPER, MAXPRM, MAX3, MXND, A(K(1)), A(K(2)),
     &            IA(K(3)), IA(K(5)), A(K(6)), A(K(7)), A(K(8)),
     &            IA(K(9)), IA(K(10)), IA(K(11)), IA(K(12)), IA(K(13)),
     &            A(K(16)), A(K(17)), IA(K(18)), IA(K(21)), IAVAIL,
     &            NAVAIL, CCW, HALFC, ERR)

C  CALCULATE A "TRIANGULAR" MAPPED MESH

            ELSE IF (TRIANG) THEN
               CALL TMSCHM (NPER, KKK, LLL, NNN, ML, MS, NSPR(L),
     &            ISLIST(IFSIDE(L)), NINT, IFLINE, NLPS, ILLIST, LINKL,
     &            LINKS, MXNPER, MAXPRM, MAX3, MXND, A(K(1)), A(K(2)),
     &            IA(K(3)), IA(K(5)), A(K(6)), A(K(7)), A(K(8)),
     &            IA(K(9)), IA(K(10)), IA(K(11)), IA(K(12)), IA(K(13)),
     &            A(K(16)), A(K(17)), IA(K(18)), IA(K(21)), IAVAIL,
     &            NAVAIL, CCW, ERR)

C  CALCULATE A "PENTAGON" MAPPED MESH

            ELSE IF (PENTAG) THEN
               CALL UMSCHM (IA, NPER, KKK, LLL, NNN, ML, MS, NSPR(L),
     &            ISLIST(IFSIDE(L)), NINT, IFLINE, NLPS, ILLIST, LINKL,
     &            LINKS, MXNPER, MAXPRM, MAX3, MXND, A(K(1)), A(K(2)),
     &            IA(K(3)), IA(K(5)), A(K(6)), A(K(7)), A(K(8)),
     &            IA(K(9)), IA(K(10)), IA(K(11)), IA(K(12)), IA(K(13)),
     &            A(K(16)), A(K(17)), IA(K(18)), IA(K(21)), IAVAIL,
     &            NAVAIL, CCW, ERR)

C  USE THE PAVING TECHNIQUE TO FILL THE INITIAL REGION

            ELSE IF (FILL) THEN
               CALL PMSCHM (NPER, NPRM, MXND, MLN, MP, ML, MS, MR, NL,
     &            MXNL, MXNPER, MAXPRM, MAXNB, MAXNBC, MAXSBC, KNBC,
     &            KSBC, KNUM, IPOINT, COOR, IPBOUN, ILINE, LTYPE, NINT,
     &            FACTOR, LCON, ILBOUN, ISBOUN, ISIDE, NLPS, IFLINE,
     &            ILLIST, ISLIST, IREGN, NPPF, IFPB, LISTPB, NLPF, IFLB,
     &            LISTLB, NSPF, IFSB, LISTSB, LINKP, LINKL, LINKS,
     &            LINKR, LINKPB, LINKLB, LINKSB, NSPR, IFSIDE, RSIZE,
     &            IFHOLE, NHPR, IHLIST, A(K(1)), A(K(2)), IA(K(3)),
     &            IA(K(4)), A(K(7)), A(K(8)), A(K(30)), IA(K(9)),
     &            IA(K(10)), IA(K(11)), IA(K(12)), IA(K(13)), IA(K(14)),
     &            IA(K(20)), A(K(22)), A(K(23)), IA(K(24)), IA(K(25)),
     &            IA(K(26)), IA(K(27)), IA(K(28)), IA(K(29)), KKK, NNN,
     &            LLL, IAVAIL, NAVAIL, DEV1, ABS(IREGN(L)), L, BATCH,
     &            NOROOM, ERR, AMESUR, XNOLD, YNOLD, NXKOLD, MMPOLD,
     &            LINKEG, LISTEG, BMESUR, MLINK, NPROLD, NPNOLD, NPEOLD,
     &            NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS,
     &            SIZMIN, EMAX, EMIN)
               IF (NOROOM) THEN
                  MXND = INT(MXND*1.5 + 1)
                  MAXNBC = MAX0 (MAXNBC, KNBC)
                  MAXSBC = MAX0 (MAXSBC, KSBC)
                  KSBC = KKSBC
                  CALL MDLONG ('XN', K(7), MXND)
                  CALL MDLONG ('YN', K(8), MXND)
                  CALL MDLONG ('ZN', K(30), MXND)
                  CALL MDLONG ('NUID', K(9), MXND)
                  CALL MDLONG ('LXK', K(10), MXND*4)
                  CALL MDLONG ('KXL', K(11), MXND*6)
                  CALL MDLONG ('NXL', K(12), MXND*6)
                  CALL MDLONG ('LXN', K(13), MXND*4)
                  CALL MDLONG ('LSTNBC', K(14), MAXNBC)
                  CALL MDLONG ('LSTSBC', K(15), MAXSBC)
                  CALL MDLONG ('NXH', K(19), MXND)
                  CALL MDLONG ('INDX', K(21), MXND)
                  CALL MDLONG ('FANGLE', K(22), MXND)
                  CALL MDLONG ('BNSIZE', K(23), MXND * 2)
                  CALL MDLONG ('LNODES', K(24), MXND * MLN)
                  CALL MDSTAT (NERR, MUSED)
                  IF (NERR .GT. 0) THEN
                     CALL MDEROR (6)
                     STOP ' '
                  END IF
                  CALL MESAGE
     &               ('REDIMENSIONING NEEDED - PLEASE WAIT')
                  IF (STEP) THEN
                     CALL MESAGE
     &                  ('CURRENT PROCESSING SCHEME IS SAVED')
                  ELSE
                     CALL MESAGE
     &                  ('CURRENT SCHEME WILL BE REPEATED')
                  END IF
                  GO TO 330
               END IF

C  PROCESS A "NORMAL" REGION

            ELSE
               CALL MMSCHM (NPER, KKK, LLL, NNN, ML, MS, NSPR(L),
     &            ISLIST(IFSIDE(L)), NINT, IFLINE, NLPS, ILLIST, LINKL,
     &            LINKS, MXNPER, MAXPRM, MAX3, MXND, A(K(1)), A(K(2)),
     &            IA(K(3)), IA(K(5)), A(K(6)), A(K(7)), A(K(8)),
     &            IA(K(9)), IA(K(10)), IA(K(11)), IA(K(12)), IA(K(13)),
     &            IAVAIL, NAVAIL, CCW, REAL, SCHSTR, M1, ERR)
            END IF

C  FLAG THE REGION IF AN ERROR HAS OCCURRED

            IF (ERR) THEN
               CALL MESAGE ('ERROR IN INITIAL GRID GENERATION')
               CALL MESAGE ('** REGION PROCESSING ABORTED **')
               CALL MESAGE (' ')
               CALL PLTBEL
               CALL PLTFLU
               IREGN(L) = ABS(IREGN(L))
               DO 340 J = 1, N(9)
                  IF (IRPB(J) .EQ. IREGN(L)) THEN
                     ADDLNK = .TRUE.
                     IMINUS = -L
                     CALL LTSORT (MR, LINKR, IREGN(L), IMINUS,
     &                  ADDLNK)
                     ADDLNK = .FALSE.
                  END IF
  340          CONTINUE
               GO TO 350
            END IF

C  BEGIN FULL SCHEME CONTROL

            RECT = .NOT.(PENTAG .OR. TRIANG .OR. TRNSIT .OR. FILL)
            IF (STEP) CALL MINMAX_FQ (MXNPER, NPER, A(K(1)), A(K(2)),
     *        XMIN, XMAX, YMIN, YMAX)
            NNNOLD = 0
            KKKOLD = 0

            CALL PSCHEM (MP, ML, MS, MR, N, IPOINT, COOR, IPBOUN, ILINE,
     &         LTYPE, NINT, FACTOR, LCON, ILBOUN, ISBOUN, ISIDE, NLPS,
     &         IFLINE, ILLIST, IREGN, NSPR, IFSIDE, ISLIST, NPPF, IFPB,
     &         LISTPB, NLPF, IFLB, LISTLB, NSPF, IFSB, LISTSB, LINKP,
     &         LINKL, LINKS, LINKR, LINKSC, LINKPB, LINKLB, LINKSB,
     &         IFHOLE, NHPR, IHLIST, MAXNBC, KNBC, MAXSBC, KSBC, MXND,
     &         NNN, NNNOLD, KKK, KKKOLD, LLL, A(K(1)), A(K(2)),
     &         IA(K(3)), IA(K(4)), A(K(7)), A(K(8)), IA(K(9)),
     &         IA(K(10)), IA(K(11)), IA(K(12)), IA(K(13)), IA(K(14)),
     &         IA(K(19)), IA(K(20)), IA(K(26)), IAVAIL, NAVAIL, MXNL,
     &         MXNPER, NPER, MAXPRM, NPRM, MSC, ISCHM, SCHEME, SCHSTR,
     &         RECT, M1, INSIDE, JJHOLE, KKSBC, DEV1, EIGHT, NINE, STEP,
     &         L, NL, MCOM, CIN, IIN, RIN, KIN, ICOM, JCOM, XMIN, XMAX,
     &         YMIN, YMAX, ICODE, NOROOM, ERR, AMESUR, XNOLD, YNOLD,
     &         NXKOLD, MMPOLD, LINKEG, LISTEG, BMESUR, MLINK, NPROLD,
     &         NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX, REYMIN,
     &         REYMAX, IDIVIS, SIZMIN, EMAX, EMIN)
            IF (NOROOM) THEN
               MXND = INT(MXND*1.2 + 1)
               MAXNBC = 2*MAX0 (MAXNBC, KNBC)
               MAXSBC = 2*MAX0 (MAXSBC, KSBC)
               KSBC = KKSBC
               CALL MDLONG ('XN', K(7), MXND)
               CALL MDLONG ('YN', K(8), MXND)
               CALL MDLONG ('ZN', K(30), MXND)
               CALL MDLONG ('NUID', K(9), MXND)
               CALL MDLONG ('LXK', K(10), MXND*4)
               CALL MDLONG ('KXL', K(11), MXND*6)
               CALL MDLONG ('NXL', K(12), MXND*6)
               CALL MDLONG ('LXN', K(13), MXND*4)
               CALL MDLONG ('LSTNBC', K(14), MAXNBC)
               CALL MDLONG ('LSTSBC', K(15), MAXSBC)
               CALL MDLONG ('NXH', K(19), MXND)
               CALL MDLONG ('INDX', K(21), MXND)
               CALL MDLONG ('FANGLE', K(22), MXND)
               CALL MDLONG ('BNSIZE', K(23), MXND * 2)
               CALL MDLONG ('LNODES', K(24), MXND * MLN)
               CALL MDSTAT (NERR, MUSED)
               IF (NERR .GT. 0) THEN
                  CALL MDEROR (6)
                  STOP ' '
               END IF
               CALL MESAGE
     &            ('REDIMENSIONING NEEDED - PLEASE WAIT')
               IF (STEP) THEN
                  CALL MESAGE ('CURRENT PROCESSING SCHEME IS SAVED')
               ELSE
                  CALL MESAGE ('CURRENT SCHEME WILL BE REPEATED')
               END IF
               GO TO 330
            ELSE IF (ERR) THEN
               IF (STEP) THEN
                  CALL INTRUP ('WOULD YOU LIKE TO REPROCESS REGION',
     &               IANS, MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN)
                  IF (IANS) GO TO 330
               END IF
               CALL MESAGE ('REGION PROCESSING ABORTED')
               GO TO 350
            ELSE IF (ICODE .EQ. IEXIT) THEN
               CALL CHKKXL (MXND, IA(K(10)), IA(K(11)), LLL, ERR)
               IF (ERR) THEN
                  CALL MESAGE ('ERROR IN CHECK OF KXL ARRAY')
                  IF (STEP) THEN
                     CALL INTRUP ('WOULD YOU LIKE TO REPROCESS REGION',
     &                  IANS, MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN)
                     IF (IANS) GO TO 330
                  END IF
                  CALL MESAGE ('REGION PROCESSING ABORTED')
                  GO TO 350
               END IF
               BAR = .FALSE.
               KSBC = 0
               CALL GETSBC (MXND, MXNPER, NPER, NL, ML, MAXSBC, MAXPRM,
     &            NPRM, IA(K(3)), IA(K(4)), A(K(7)), A(K(8)), IA(K(9)),
     &            IA(K(10)), IA(K(11)), IA(K(12)), IA(K(15)), IA(K(20)),
     &            KSBC, LCON, ISBOUN, LINKL, NSPF, IFSB, LISTSB, LINKSB,
     &            LLL, BAR, ERR)
               IF (ERR) THEN
                  CALL MESAGE ('ERROR IN SORTING SIDE BOUNDARIES')
                  IF (STEP) THEN
                     CALL INTRUP ('WOULD YOU LIKE TO REPROCESS REGION',
     &                  IANS, MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN)
                     IF (IANS) GO TO 330
                  END IF
                  CALL MESAGE ('REGION PROCESSING ABORTED')
                  GO TO 350
               END IF
               CALL MKUSED (MXNL, MP, ML, IA(K(4)), IPOINT, NINT, LINKP,
     &            LINKL, LCON, NL)
               CALL SAVREG (MXND, MAXNBC, MAXSBC, A(K(7)), A(K(8)),
     &            IA(K(9)), IA(K(10)), IA(K(12)), IA(K(13)), IA(K(14)),
     &            IA(K(15)), KNBC, KSBC, NNN, KKK, ABS(IREGN(L)), IUNIT,
     &            BAR, M1)
               IREGN(L) = ABS(IREGN(L))
               NPREGN = NPREGN + 1
               NPELEM = NPELEM + KKK
               NPNODE = NPNODE + NNN
               IF (EIGHT .OR. NINE) NPNODE = NPNODE + LLL
               IF (NINE) NPNODE = NPNODE + KKK
               NPNBC = NPNBC + KNBC
               IF (THREE .OR. EIGHT .OR. NINE)
     &            NPNBC = NPNBC + KNBC
               NPSBC = NPSBC + KSBC
               IF (THREE .OR. EIGHT .OR. NINE)
     &            NPSBC = NPSBC + KSBC
               MAXKXN = MAXKXN + LLL
               IREGN(L) = ABS(IREGN(L))
               WRITE (*, 10110) IREGN(L)
            ELSE IF (ICODE .EQ. IOVER) THEN
               GO TO 330
            ELSE IF (ICODE .EQ. IQUIT) THEN
               GO TO 350
            END IF
         END IF
         END IF
  350 CONTINUE

C  END OF THIS SET OF REGIONS
C     IF STEPPING THROUGH, SEE IF ANY MORE REGIONS
C        ARE TO BE PROCESSED

  360 CONTINUE
      IF (STEP) THEN
         CALL INTRUP ('PROCESS ADDITIONAL REGIONS', IANS, MCOM,
     &      ICOM, JCOM, CIN, IIN, RIN, KIN)
         IF (IANS) GO TO 320
      END IF

C  SET UP THE LOOP FOR PROCESSING BAR SETS

  370 CONTINUE
  380 CONTINUE
      IF (STEP .AND. (N(21) .GT. 0)) THEN
         CALL MESAGE ('STEP PROCESS BAR SETS I1 THROUGH I2')
         IF (ICOM .GT. JCOM) THEN
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN, IIN,
     &         RIN)
            ICOM = 1
         END IF
         CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2, IFOUND)
         IF (IFOUND .GT. 0) THEN
            CALL CHECK (I1, I2, N(21))
         ELSE
            GO TO 450
         END IF
      ELSE
         I1 = 1
         I2 = N(21)
      END IF

C  BEGIN PROCESSING BAR SETS

      REAL = .TRUE.
      COUNT = .FALSE.
      DO 440 I = I1, I2
         CALL LTSORT (MS, LINKB, I, IPNTR, ADDLNK)

C  SEE IF THIS BAR SET IS FOR SPRINGS

         IF ((IPNTR .GT. 0) .AND. (IBARST(IPNTR) .LT. 0) .AND.
     &      (JMAT(IPNTR) .LT. 0)) THEN
            L = IPNTR
            WRITE (*, 10130) ABS(IBARST(L))

            CALL SPRING (MP, ML, MS, MXNPER, MXND, MAXNBC, MAXSBC, L,
     &         IPOINT, COOR, IPBOUN, LINKP, ILINE, LTYPE, NINT, FACTOR,
     &         LCON, ILBOUN, ISBOUN, LINKL, NLPB, JFLINE, JLLIST,
     &         LINKPB, NPPF, IFPB, LISTPB, LINKLB, NLPF, IFLB, LISTLB,
     &         LINKSB, NSPF, IFSB, LISTSB, LSTNBC, A(K(1)), A(K(2)),
     &         IA(K(3)), A(K(7)), A(K(8)), IA(K(9)), IA(K(10)), NNN,
     &         KKK, LLL, KNBC, KSBC, ERR, ADDLNK, COUNT, NOROOM,
     &         AMESUR, XNOLD, YNOLD, NXKOLD, MMPOLD, LINKEG, LISTEG,
     &         BMESUR, MLINK, NPROLD, NPNOLD, NPEOLD, NNXK, REMESH,
     &         REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX,
     &         EMIN, GRAPH)
            IF (ERR) THEN
               CALL MESAGE ('ERROR IN 2-NODE SPRING ELEMENT '//
     &            'GENERATION')
               CALL MESAGE ('BAR SET PROCESSING ABORTED')
               CALL MESAGE (' ')
               CALL PLTBEL
               CALL PLTFLU
               IBARST(L) = ABS(IBARST(L))
               DO 390 M = 1, N(9)
                  IF (ABS(IRPB(M)) .EQ. IBARST(L)) THEN
                     CALL LTSORT (MS, LINKB, IRPB(M), IPNTR, ADDLNK)
                     ADDLNK = .TRUE.
                     IMINUS = -IPNTR
                     CALL LTSORT (MS, LINKB, IRPB(M), IMINUS, ADDLNK)
                     ADDLNK = .FALSE.
                  END IF
  390          CONTINUE
               GO TO 430
            END IF

C  PROCESS A REGULAR BARSET

         ELSE IF ((IPNTR .GT. 0) .AND. (IBARST(IPNTR) .LT. 0)) THEN
            L = IPNTR
            WRITE (*, 10120) ABS(IBARST(L))
            REAL = .TRUE.
            TEST = .FALSE.
            KKK = 0
            NNN = 0
            KNBC = 0
            KSBC = 0
            LLL = 1

C  LOOP THROUGH ALL THE LINES IN THE BAR SETS

            DO 410 J = JFLINE(L), JFLINE(L) + NLPB(L) - 1
               CALL LTSORT (ML, LINKL, JLLIST(J), KK, ADDLNK)
               CALL LTSORT (MP, LINKP, LCON(1, KK), IP1, ADDLNK)
               CALL LTSORT (MP, LINKP, LCON(2, KK), IP2, ADDLNK)
               IF (LCON(3, KK) .GT. 0) THEN
                  CALL LTSORT (MP, LINKP, LCON(3, KK), IP3, ADDLNK)
               ELSE IF (LCON(3, KK) .LT. 0) THEN
                  CALL LTSORT (MP, LINKP, ABS(LCON(3, KK)), IP3,
     &               ADDLNK)
                  IP3 = -IP3
               ELSE
                  IP3 = 0
               END IF

C  CALCULATE NODES IN THE BAR SET LINE

               CALL PLINE (MP, ML, MXNPER, MAXNBC, MAXSBC, IPOINT,
     &            COOR, LINKP, ILINE(KK), LTYPE(KK), NINT(KK),
     &            FACTOR(KK), IP1, IP2, IP3, A(K(1)), A(K(2)), IA(K(3)),
     &            IPBOUN(IP1), IPBOUN(IP2), ILBOUN(KK), ISBOUN(KK),
     &            LINKPB, NPPF, IFPB, LISTPB, LINKLB, NLPF, IFLB,
     &            LISTLB, LINKSB, NSPF, IFSB, LISTSB, IA(K(14)), KNBC,
     &            KSBC, ERR, TEST, REAL, COUNT, NOROOM, AMESUR, XNOLD,
     &            YNOLD, NXKOLD, MMPOLD, LINKEG, LISTEG, BMESUR, MLINK,
     &            NPROLD, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN, REXMAX,
     &            REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN, GRAPH,
     &            DXMAX)
               IF (ERR) THEN
                  CALL MESAGE ('ERROR IN 2-NODE ELEMENT GENERATION')
                  CALL MESAGE ('BAR SET PROCESSING ABORTED')
                  CALL MESAGE (' ')
                  CALL PLTBEL
                  CALL PLTFLU
                  IBARST(L) = ABS(IBARST(L))
                  DO 400 M = 1, N(9)
                     IF (ABS(IRPB(M)) .EQ. IBARST(L)) THEN
                        CALL LTSORT (MS, LINKB, IRPB(M), IPNTR, ADDLNK)
                        ADDLNK = .TRUE.
                        IMINUS = -IPNTR
                        CALL LTSORT (MS, LINKB, IRPB(M), IMINUS, ADDLNK)
                        ADDLNK = .FALSE.
                     END IF
  400             CONTINUE
                  GO TO 430
               END IF

C  ADD THESE NODES AND ELEMENTS TO THE CURRENT LIST

               NNN0 = NNN + 1
               NNN = NNN + ABS(NINT(KK)) + 1
               IF (JCENT(L) .GT. 0) THEN
                  CALL LTSORT (MP, LINKP, JCENT(L), IP3, ADDLNK)
               ELSE
                  IP3 = 0
               END IF
               CALL MAK2EL (MP, MXNPER, MXND, NNN0, NNN, KKK, A(K(1)),
     &            A(K(2)), IA(K(3)), A(K(7)), A(K(8)), IA(K(9)),
     &            IA(K(10)), COOR, IP3)

C  MARK THESE POINTS AND THE LINE AS BEING USED

               NINT(KK) = -ABS(NINT(KK))
               IPOINT(IP1) = -ABS(IPOINT(IP1))
               IPOINT(IP2) = -ABS(IPOINT(IP2))
  410       CONTINUE
         ENDIF

C  WRITE OUT THE BAR SET ELEMENTS AND BOUNDARY CONDITIONS

         IF ((IPNTR .GT. 0) .AND. (IBARST(IPNTR) .LT. 0)) THEN
            BAR = .TRUE.
            KSBC = 0
            CALL GETSBC (MXND, NNN, NNN, NLPB(L), ML, MAXSBC, 1,
     &         1, IA(K(3)), JLLIST(JFLINE(L)), A(K(7)), A(K(8)),
     &         IA(K(9)), IA(K(10)), IA(K(11)), IA(K(12)), IA(K(15)),
     &         IA(K(20)), KSBC, LCON, ISBOUN, LINKL, NSPF, IFSB,
     &         LISTSB, LINKSB, KKK, BAR, ERR)
            IF (ERR) THEN
               CALL MESAGE ('ERROR IN SORTING SIDE BOUNDARIES')
               CALL MESAGE ('BAR SET PROCESSING ABORTED')
               CALL MESAGE (' ')
               IBARST(L) = ABS(IBARST(L))
               DO 420 M = 1, N(9)
                  IF (ABS(IRPB(M)) .EQ. IBARST(L)) THEN
                     CALL LTSORT (MS, LINKB, IRPB(M), IPNTR, ADDLNK)
                     ADDLNK = .TRUE.
                     IMINUS = -IPNTR
                     CALL LTSORT (MS, LINKB, IRPB(M), IMINUS, ADDLNK)
                     ADDLNK = .FALSE.
                  END IF
  420          CONTINUE
               GO TO 430
            END IF
            CALL SAVREG (MXND, MAXNBC, MAXSBC, A(K(7)), A(K(8)),
     &         IA(K(9)), IA(K(10)), IA(K(12)), IA(K(13)), IA(K(14)),
     &         IA(K(15)), KNBC, KSBC, NNN, KKK, ABS(IBARST(L)), IUNIT,
     &         BAR, M1)
            NPREGN = NPREGN + 1
            NPELEM = NPELEM + KKK
            NPNODE = NPNODE + NNN
            IF (THREE) NPNODE = NPNODE + MAX0 (NNN, KKK + 1)
            NPNBC = NPNBC + KNBC
            IF (THREE) NPNBC = NPNBC + KNBC
            NPSBC = NPSBC + KSBC
            IF (THREE) NPSBC = NPSBC + KSBC
            MAXKXN = MAXKXN + KKK + NLPB(L)
            IBARST(L) = ABS(IBARST(L))
            WRITE (*, 10140) IBARST(L)
         END IF

C  END OF THIS BAR SET

  430    CONTINUE
  440 CONTINUE

C  END OF THIS GROUP OF BAR SETS
C  IF STEPPING THROUGH, SEE IF ANY MORE BAR SETS ARE TO BE PROCESSED

      IF (STEP .AND. (N(21) .GT. 0)) THEN
         IF ((ICOM .LE. JCOM) .AND. ((CIN(ICOM)(1:1) .EQ. 'Y') .OR.
     &      (CIN(ICOM)(1:1) .EQ. 'y'))) THEN
            IANS = .TRUE.
            ICOM = ICOM + 1
         ELSE IF ((ICOM .LE. JCOM) .AND. ((CIN(ICOM)(1:1) .EQ. 'N') .OR.
     &      (CIN(ICOM)(1:1) .EQ. 'n'))) THEN
            IANS = .TRUE.
            ICOM = ICOM + 1
         ELSE
            IF ((ICOM .LE. JCOM) .AND. ((CIN(ICOM)(1:1) .EQ. 'Y') .OR.
     &         (CIN(ICOM)(1:1) .EQ. 'y'))) THEN
               IANS = .TRUE.
               ICOM = ICOM + 1
            ELSE IF ((ICOM .LE. JCOM) .AND. ((CIN(ICOM)(1:1) .EQ. 'N')
     &         .OR. (CIN(ICOM)(1:1) .EQ. 'n'))) THEN
               IANS = .FALSE.
               ICOM = ICOM + 1
            ELSE
               CALL MESAGE (' ')
               CALL INTRUP ('PROCESS ADDITIONAL BAR SETS', IANS, MCOM,
     &            ICOM, JCOM, CIN, IIN, RIN, KIN)
            END IF
         END IF
         IF (IANS) GO TO 380
      END IF

C  RESTORE THE DATA BASE TO ITS INITIAL CONDITION

  450 CONTINUE
      DO 460 I = 1, N(1)
         IPOINT(I) = ABS(IPOINT(I))
  460 CONTINUE
      DO 470 I = 1, N(2)
         NINT(I) = ABS(NINT(I))
  470 CONTINUE
      DO 480 I = 1, N(5)
         IBARST(I) = ABS(IBARST(I))
  480 CONTINUE
      DO 490 I = 1, N(7)
         IREGN(I) = ABS(IREGN(I))
  490 CONTINUE
      DO 500 I = 1, N(5)
         CALL LTSORT (MS, LINKB, IBARST(I), IPNTR, ADDLNK)
         ADDLNK = .TRUE.
         IPLUS = ABS(IPNTR)
         CALL LTSORT (MS, LINKB, IBARST(I), IPLUS, ADDLNK)
         ADDLNK = .FALSE.
  500 CONTINUE
      DO 510 I = 1, N(7)
         CALL LTSORT (MR, LINKR, IREGN(I), IPNTR, ADDLNK)
         ADDLNK = .TRUE.
         IPLUS = ABS(IPNTR)
         CALL LTSORT (MR, LINKR, IREGN(I), IPLUS, ADDLNK)
         ADDLNK = .FALSE.
  510 CONTINUE

      CALL MDDEL ('X')
      CALL MDDEL ('Y')
      CALL MDDEL ('NID')
      CALL MDDEL ('LISTL')
      CALL MDDEL ('NNPS')
      CALL MDDEL ('ANGLE')
      CALL MDDEL ('XN')
      CALL MDDEL ('YN')
      CALL MDDEL ('NUID')
      CALL MDDEL ('LXK')
      CALL MDDEL ('KXL')
      CALL MDDEL ('NXL')
      CALL MDDEL ('LXN')
      CALL MDDEL ('LSTNBC')
      CALL MDDEL ('LSTSBC')
      CALL MDDEL ('XSUB')
      CALL MDDEL ('YSUB')
      CALL MDDEL ('NIDSUB')
      CALL MDDEL ('NXH')
      CALL MDDEL ('NPERIM')
      CALL MDDEL ('INDX')
      CALL MDDEL ('FANGLE')
      CALL MDDEL ('BNSIZE')
      CALL MDDEL ('LNODES')
      CALL MDDEL ('PRLINK')
      CALL MDDEL ('MARKED')
      CALL MDDEL ('IPTPER')
      CALL MDDEL ('NUMPER')
      CALL MDDEL ('LPERIM')
      CALL MDDEL ('ZN')
      CALL MDSTAT (NERR, MUSED)
      IF (NERR .GT. 0) THEN
         CALL MDEROR (6)
         STOP ' '
      END IF
      RETURN

10000 FORMAT (' INITIAL CHECK BEGUN FOR REGION:', I5)
10010 FORMAT (' INITIAL CHECK BEGUN FOR GROUP:', I5)
10020 FORMAT (' ...INITIAL CHECK BEGUN FOR REGION:', I5)
10030 FORMAT (' ** ERROR - REGION', I5, ' IN BODY LIST IS INVALID **')
10040 FORMAT (' INITIAL CHECK BEGUN FOR BAR SET:', I5)
10050 FORMAT (' ** ERROR - LINE GENERATION ERRORS FOR BAR  SET', I5,
     &   ' **')
10060 FORMAT (' ** ERROR - BAR SET', I5, ' IN BODY LIST IS INVALID **')
10070 FORMAT (' NOW PROCESSING GROUP: ', I5)
10080 FORMAT (' ...NOW PROCESSING REGION:', I5)
10090 FORMAT (' NOW PROCESSING REGION:', I5)
10100 FORMAT (' GROUP', I5, ' SUCCESSFULLY COMPLETED AND SAVED')
10110 FORMAT (' REGION', I5, ' SUCCESSFULLY COMPLETED AND SAVED')
10120 FORMAT (' NOW PROCESSING BAR SET:', I5)
10130 FORMAT (' NOW PROCESSING SPRING BAR SET:', I5)
10140 FORMAT (' BAR SET', I5, ' SUCCESSFULLY COMPLETED AND SAVED')
      END
