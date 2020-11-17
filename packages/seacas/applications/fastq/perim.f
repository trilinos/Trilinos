C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE PERIM (MP, ML, MS, NS, MAXNL, MAXNP, MAXNBC, MAXSBC,
     &   KNBC, KSBC, KNUM, IPOINT, COOR, IPBOUN, ILINE, LTYPE, NINT,
     &   FACTOR, LCON, ILBOUN, ISBOUN, ISIDE, NLPS, IFLINE, ILLIST,
     &   ISLIST, NPPF, IFPB, LISTPB, NLPF, IFLB, LISTLB, NSPF, IFSB,
     &   LISTSB, LINKP, LINKL, LINKS, LINKPB, LINKLB, LINKSB, X, Y, NID,
     &   N, LISTL, NL, LSTNBC, MARKED, EVEN, REAL, ERR, CCW, COUNT,
     &   NOROOM, AMESUR, XNOLD, YNOLD, NXKOLD, MMPOLD, LINKEG, LISTEG,
     &   BMESUR, MLINK, NPROLD, NPNOLD, NPEOLD, NNXK, REMESH, REXMIN,
     &   REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN)
C***********************************************************************

C  SUBROUTINE PERIM = GENERATES THE PERIMETER OF A REGION

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     QMESH = GENERATES THE QUADRILATERAL MESH

C***********************************************************************

C  VARIABLES USED:
C     X   = THE X VALUES OF THE PERIMETER LIST
C     Y   = THE Y VALUES OF THE PERIMETER LIST
C     NID = AN ARRAY OF UNIQUE NODE IDENTIFIERS WITHIN THE BODY.
C           IF 00000YYYYY,  YYYYY IS AN INDEX INTO THE POINT TABLE.
C           IF 1XXXXYYYYY,  XXXX IS AN INDEX INTO THE LINE TABLE.
C     N   = THE NUMBER OF NODES ON THE PERIMETER
C     ERR = .TRUE. IF ERRORS WERE ENCOUNTERED

C***********************************************************************

      DIMENSION IPOINT (MP), COOR (2, MP), IPBOUN (MP)
      DIMENSION ILINE (ML), NINT (ML), LTYPE (ML)
      DIMENSION FACTOR (ML), LCON (3, ML)
      DIMENSION ILBOUN (ML), ISBOUN (ML)
      DIMENSION ISIDE (MS), NLPS (MS), IFLINE (MS)
      DIMENSION ILLIST (MS*3), ISLIST (NS)
      DIMENSION LINKP (2, MP), LINKL (2, ML), LINKS (2, MS)
      DIMENSION LISTL (MAXNL), X (MAXNP), Y (MAXNP), NID (MAXNP)
      DIMENSION LINKPB (2, MP), NPPF (MP), IFPB (MP), LISTPB (2, MP)
      DIMENSION LINKLB (2, ML), NLPF (ML), IFLB (ML), LISTLB (2, ML)
      DIMENSION LINKSB (2, ML), NSPF (ML), IFSB (ML), LISTSB (2, ML)
      DIMENSION LSTNBC (MAXNBC), MARKED (3, MAXNL)

      DIMENSION AMESUR(NPEOLD), XNOLD(NPNOLD), YNOLD(NPNOLD)
      DIMENSION NXKOLD(NNXK, NPEOLD), MMPOLD(3, NPROLD)
      DIMENSION LINKEG(2, MLINK), LISTEG(4 * NPEOLD), BMESUR(NPNOLD)

      LOGICAL ERR, REAL, EVEN, CCW, TEST, ADDLNK, ITRIED, INDETR, NOROOM
      LOGICAL COUNT, REMESH, GRAPH

      N = 0
      ADDLNK = .FALSE.
      GRAPH = .FALSE.
      ITRIED = .FALSE.
      IF (GRAPH) THEN
         XDIMD = 1.
         YDIMD = .75
         CALL PLTBGN
         XDIMR = REXMAX - REXMIN
         YDIMR = REYMAX - REYMIN
         CALL MPVIEW (0., XDIMD, 0., YDIMD)
         XRAT = XDIMR/XDIMD
         YRAT = YDIMR/YDIMD
         IF (XRAT.LT.YRAT) THEN
            XDIMR = XDIMD * YRAT
            XX1 =  (REXMIN + REXMAX - XDIMR) * .5
            XX2 =  (REXMIN + REXMAX + XDIMR) * .5
            XDIMR = XX2 - XX1
            YY1 = REYMIN
            YY2 = REYMAX
         ELSE
            YDIMR = YDIMD * XRAT
            YY1 =  (REYMIN + REYMAX - YDIMR) * .5
            YY2 =  (REYMIN + REYMAX + YDIMR) * .5
            YDIMR = YY2 - YY1
            XX1 = REXMIN
            XX2 = REXMAX
         ENDIF
         XX1 = XX1 - (XDIMR * .1)
         XX2 = XX2 + (XDIMR * .1)
         YY1 = YY1 - (YDIMR * .1)
         YY2 = YY2 + (YDIMR * .1)
         CALL MPORT2 (XX1, XX2, YY1, YY2)
         CALL PLTFRM (0)
      ENDIF

C  GET LIST OF LINES

      CALL LLIST (MS, ML, MAXNL, NS, NL, KNUM, LISTL, ILINE, ISIDE,
     &   NLPS, IFLINE, ILLIST, LCON, ISLIST, LINKS, LINKL, ERR)
      IF (ERR)RETURN
      ERR = .TRUE.

C  DEFINE VALUE OF KP,  THE BEGINNING CONNECTIVITY POINT

      IF (NL .LT. 2) THEN
         CALL LTSORT (ML, LINKL, LISTL (1), ILI, ADDLNK)
         KP = LCON (1, ILI)
      ELSE
         CALL LTSORT (ML, LINKL, LISTL (1), IL1, ADDLNK)
         CALL LTSORT (ML, LINKL, LISTL (2), IL2, ADDLNK)
         J1 = LCON (1, IL1)
         J2 = LCON (2, IL1)
         K1 = LCON (1, IL2)
         K2 = LCON (2, IL2)
         KP = 0
         IF ((J1 .EQ. K1).OR. (J1 .EQ. K2))KP = J2
         IF ((J2 .EQ. K1).OR. (J2 .EQ. K2))KP = J1
         ILI = IL2
         IF (KP .EQ. 0) THEN
            WRITE (*, 10000)ILINE (IL2)
            RETURN
         ENDIF
      ENDIF

C  IF PERIMETER HAS ODD NUMBER OF POINTS,
C  DO A TRIAL PERIMETER GENERATION TO SEE WHERE TO INSERT A NODE.

      IF (EVEN) THEN
         NUMINT = 0
         DO 100 IL = 1, NL
            CALL LTSORT (ML, LINKL, LISTL (IL), ILI, ADDLNK)
            NUMINT = NUMINT+IABS (NINT (ILI))
  100    CONTINUE
         IF ((MOD (NUMINT, 2) .EQ. 1) .OR. REMESH) THEN
            DX = 0.0
            IX = 0
            DO 140 IL = 1, NL
               CALL LTSORT (ML, LINKL, LISTL (IL), ILI, ADDLNK)

C  SKIP PREVIOUSLY USED LINES
C  SKIP LINES USED TWICE IN THIS BOUNDARY  (RECALL ANNULUS)

               IF (NINT (ILI) .GT. 0) THEN
                  DO 110 I = 1, NL
                     CALL LTSORT (ML, LINKL, LISTL (I), IPNTR, ADDLNK)
                     IF ((I.NE.IL).AND. (IPNTR .EQ. ILI))GOTO 130
  110             CONTINUE
                  CALL LTSORT (MP, LINKP, LCON (1, ILI), IP1, ADDLNK)
                  CALL LTSORT (MP, LINKP, LCON (2, ILI), IP2, ADDLNK)
                  IF (LCON (3, ILI) .GT. 0) THEN
                     CALL LTSORT (MP, LINKP, LCON (3, ILI), IP3, ADDLNK)
                  ELSEIF (LCON (3, ILI) .LT. 0) THEN
                     CALL LTSORT (MP, LINKP, IABS (LCON (3, ILI)), IP3,
     &                  ADDLNK)
                     IP3 = -IP3
                  ELSE
                     IP3 = 0
                  ENDIF
                  TEST = .TRUE.

                  CALL PLINE (MP, ML, MAXNP, MAXNBC, MAXSBC, IPOINT,
     &               COOR, LINKP, ILINE (ILI), LTYPE (ILI), NINT (ILI),
     &               FACTOR (ILI), IP1, IP2, IP3, X, Y, NID,
     &               IPBOUN (IP1), IPBOUN (IP2), ILBOUN (ILI),
     &               ISBOUN (ILI), LINKPB, NPPF, IFPB, LISTPB, LINKLB,
     &               NLPF, IFLB, LISTLB, LINKSB, NSPF, IFSB, LISTSB,
     &               LSTNBC, KNBC, KSBC, ERR, TEST, REAL, COUNT, NOROOM,
     &               AMESUR, XNOLD, YNOLD, NXKOLD, MMPOLD, LINKEG,
     &               LISTEG, BMESUR, MLINK, NPROLD, NPNOLD, NPEOLD,
     &               NNXK, REMESH, REXMIN, REXMAX, REYMIN, REYMAX,
     &               IDIVIS, SIZMIN, EMAX, EMIN, GRAPH, DXMAX)
                  IF (ERR) RETURN
                  IF (REMESH) THEN
                     IF (IL .EQ. 1) THEN
                        IX = ILI
                        DX = DXMAX
                     ELSE
                        IF (DXMAX .GT. DX) THEN
                           IX = ILI
                           DX = DXMAX
                        ENDIF
                     ENDIF
                  ELSE
                     DO 120 I = 2, IABS (NINT (ILI))+1
                        D = SQRT ((X (I)-X (I-1))**2 +
     &                     (Y (I)-Y (I-1))**2)
                        IF (D .GT. DX) THEN
                           DX = D
                           IX = ILI
                        ENDIF
  120                CONTINUE
                  ENDIF
               ENDIF
  130          CONTINUE
  140       CONTINUE
            IF (IX .EQ. 0) THEN
               ERR = .TRUE.
               WRITE (*, 10010)KNUM
               RETURN
            ENDIF

C  RECALCULATE THE NUMBER OF INTERVALS IF REMESHING

            IF (REMESH) THEN
               NUMINT = 0
               DO 150 IL = 1, NL
                  CALL LTSORT (ML, LINKL, LISTL (IL), ILI, ADDLNK)
                  NUMINT = NUMINT+IABS (NINT (ILI))
  150          CONTINUE
               IF (MOD (NUMINT, 2) .EQ. 1) THEN
                  NINT (IX) = NINT (IX) + 1
                  WRITE (*, 10020)ILINE (IX)
               ENDIF
            ELSE
               NINT (IX) = NINT (IX)+1
               WRITE (*, 10020)ILINE (IX)
            ENDIF
         ENDIF
      ENDIF

C  NOW LOOP THROUGH THE LINES TO GENERATE THE PERIMETER

      IF (GRAPH) THEN
         CALL PLTBGN
         CALL MPVIEW (0., XDIMD, 0., YDIMD)
         CALL MPORT2 (XX1, XX2, YY1, YY2)
         CALL PLTFRM (0)
      ENDIF
  160 CONTINUE
      J1 = KP
      N = 0
      DO 170 IL = 1, NL
         CALL LTSORT (ML, LINKL, LISTL (IL), ILI, ADDLNK)
         K1 = LCON (1, ILI)
         K2 = LCON (2, ILI)
         IF ((K1.NE.KP).AND. (K2.NE.KP)) THEN
            WRITE (*, 10000)ILINE (ILI)
            RETURN
         ENDIF
         IF (N .GT. 0)NIDSAV = NID (N+1)
         CALL LTSORT (MP, LINKP, LCON (1, ILI), IP1, ADDLNK)
         CALL LTSORT (MP, LINKP, LCON (2, ILI), IP2, ADDLNK)
         IF (LCON (3, ILI) .GT. 0) THEN
            CALL LTSORT (MP, LINKP, LCON (3, ILI), IP3, ADDLNK)
         ELSEIF (LCON (3, ILI) .LT. 0) THEN
            CALL LTSORT (MP, LINKP, IABS (LCON (3, ILI)), IP3, ADDLNK)
            IP3 = -IP3
         ELSE
            IP3 = 0
         ENDIF
         IMAXNP = MAXNP-N
         TEST = .FALSE.

         if (imaxnp .lt. nint(ili)) then
           stop 'ERROR: Intervals larger than space'
         end if

         CALL PLINE (MP, ML, IMAXNP, MAXNBC, MAXSBC, IPOINT, COOR,
     &      LINKP, ILINE (ILI), LTYPE (ILI), NINT (ILI), FACTOR (ILI),
     &      IP1, IP2, IP3, X (N+1), Y (N+1), NID (N+1), IPBOUN (IP1),
     &      IPBOUN (IP2), ILBOUN (ILI), ISBOUN (ILI), LINKPB, NPPF,
     &      IFPB, LISTPB, LINKLB, NLPF, IFLB, LISTLB, LINKSB, NSPF,
     &      IFSB, LISTSB, LSTNBC, KNBC, KSBC, ERR, TEST, REAL, COUNT,
     &      NOROOM, AMESUR, XNOLD, YNOLD, NXKOLD, MMPOLD, LINKEG,
     &      LISTEG, BMESUR, MLINK, NPROLD, NPNOLD, NPEOLD, NNXK, REMESH,
     &      REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX, EMIN,
     &      GRAPH, DXMAX)
         IF (ERR)RETURN
         ERR = .TRUE.

C  REVERSE LINE IF NECESSARY

         IF (K1.NE.KP) THEN
            CALL REVERS (X (N+1), IABS (NINT (ILI))+1)
            CALL REVERS (Y (N+1), IABS (NINT (ILI))+1)
            CALL IREVER (NID (N+1), IABS (NINT (ILI))+1)
         ENDIF
         IF (N .GT. 0)NID (N+1) = NIDSAV

C  FINISH UP WITH THIS LINE
C  KP IS THE POINT ON THE FAR END OF THE LINE
C  DON'T INCLUDE THE LAST POINT ON THIS LINE IN THE LIST

         KP =  (K1+K2)-KP
         N = N+IABS (NINT (ILI))

C  MARK ALL THE LINES AS USED AND
C  REMEMBER WHICH HAVE JUST BEEN MARKED

         IF (IPOINT (IP1) .GT. 0) THEN
            MARKED (1, IL) = IP1
            IPOINT (IP1) = - IABS (IPOINT (IP1))
         ELSE
            MARKED (1, IL) = 0
         ENDIF
         IF (IPOINT (IP2) .GT. 0) THEN
            MARKED (2, IL) = IP2
            IPOINT (IP2) = - IABS (IPOINT (IP2))
         ELSE
            MARKED (2, IL) = 0
         ENDIF
         IF (NINT (ILI) .GT. 0) THEN
            MARKED (3, IL) = ILI
            NINT (ILI) = - IABS (NINT (ILI))
         ELSE
            MARKED (3, IL) = 0
         ENDIF
  170 CONTINUE

C  LINES ARE EXHAUSTED  -  CHECK FOR CIRCULARITY

      CALL LTSORT (ML, LINKL, LISTL (1), ILI, ADDLNK)
      IF (KP.NE.J1) THEN
         WRITE (*, 10000)ILINE (ILI)
         RETURN
      ENDIF

C  RESET THE JUST MARKED LINES

      DO 180 I = 1, NL
         IF (MARKED (1, I) .GT. 0)
     &      IPOINT (MARKED (1, I)) = IABS (IPOINT (MARKED (1, I)))
         IF (MARKED (2, I) .GT. 0)
     &      IPOINT (MARKED (2, I)) = IABS (IPOINT (MARKED (2, I)))
         IF (MARKED (3, I) .GT. 0)
     &      NINT (MARKED (3, I)) = IABS (NINT (MARKED (3, I)))
  180 CONTINUE

C  PERIMETER COMPLETED
C  INSURE PROPER ORIENTATION  (COUNTER-CLOCKWISE)

      CALL CCLOCK (X, Y, N, CCW, ERR, INDETR)

C  THE LINE ORIENTATION MAY BE BAD - TRY A SIMPLE FIX

      IF ((INDETR).AND. (.NOT.ITRIED)) THEN
         DO 190 IL = 1, NL
            CALL LTSORT (ML, LINKL, LISTL (IL), ILI, ADDLNK)
            IF (LCON (1, ILI) .EQ. LCON (2, ILI)) THEN
               LCON (3, ILI) = -LCON (3, ILI)
               ITRIED = .TRUE.
               GOTO 160
            ENDIF
  190    CONTINUE
      ELSEIF ((INDETR).AND. (ITRIED)) THEN
         CALL MESAGE ('ORIENTATION OF PERIMETER CANNOT BE DETERMINED')
      ENDIF
      IF (ERR)RETURN
      ERR = .TRUE.
      IF (.NOT.CCW) THEN
         CALL REVERS (X (2), N-1)
         CALL REVERS (Y (2), N-1)
         CALL IREVER (NID (2), N-1)
      ENDIF

C  EXIT

      ERR = .FALSE.
      RETURN

10000 FORMAT (' LINE', I5, ' DOES NOT CONNECT TO THE', /,
     +   ' PREVIOUS SECTION OF THE PERIMETER')
10010 FORMAT (' IN REGION', I5, ' NO LINE IS ALTERABLE TO ENFORCE', /,
     +   ' AN EVEN NUMBER OF PERIMETER POINTS')
10020 FORMAT (' NO. OF INTERVALS ON LINE', I5, ' WAS INCREASED BY 1')

      END
