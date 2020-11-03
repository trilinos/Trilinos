C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE RENUM (NPNODE, NPELEM, MXNFLG, MXSFLG, NPNBC, NPSBC,
     &   NPWTS, NPREGN, MP, ML, MS, MR, MSC, MAXKXN, NNUID, NNXK,
     &   MXLPS, IUNIT, NNN, KKK, KCRD, NL, NPBF, NLBF, NSBF, IPART,
     &   LSTNBC, LSTSBC, NNFLG, NNPTR, NNLEN, NSFLG, NSPTR, NSLEN,
     &   NVPTR, NVLEN, NSIDEN, NUID, XN, YN, NXK, MAT, KXN, LIST, LA,
     &   LB, MATMAP, LISTN, WTNODE, WTSIDE, WTHOLD, IHERE, ILIST, XLIST,
     &   NUMMAT, NBCNOD, NNLIST, NBCSID, NSLIST, NVLIST, COOR, ILINE,
     &   LTYPE, LCON, ISIDE, NLPS, IFLINE, ILLIST, LINKP, LINKL, LINKS,
     &   LINKR, IMAT, LINKB, JMAT, IPBF, NPPF, IFPB, LISTPB, IWTPBF,
     &   ILBF, NLPF, IFLB, LISTLB, IWTLBF, ISBF, NSPF, IFSB, LISTSB,
     &   IWTSBF, LINKPB, LINKLB, LINKSB, NUMBER, THREE, EIGHT, NINE,
     &   OPTIM, ISBARS)
C***********************************************************************

C  SUBROUTINE RENUM = NUMBERS QMESH OUTPUT, AND RENUMBERS AS NEEDED FOR
C                     OPTIMIZATION

C***********************************************************************

C    THE REFERENCE DOCUMENTS FOR THIS CODE ARE SLA-73-1088, JULY 1974,
C    AND SLA-74-0239, JULY 1974

C***********************************************************************

      DIMENSION NLIST(20), IPART(3, NPREGN)
      DIMENSION COOR(2, MP), ILINE(ML), LTYPE(ML), LCON(3, ML)
      DIMENSION ISIDE(MS), NLPS(MS), IFLINE(MS), ILLIST(MS*3)
      DIMENSION LINKP(2, MP), LINKL(2, ML), LINKS(2, MS), LINKB(2, MS)
      DIMENSION LINKR(2, MR)
      DIMENSION IMAT(MR), JMAT(MS)
      DIMENSION IPBF(MP), NPPF(MP), IFPB(MP), LISTPB(2, MP)
      DIMENSION ILBF(ML), NLPF(ML), IFLB(ML), LISTLB(2, ML)
      DIMENSION ISBF(ML), NSPF(ML), IFSB(ML), LISTSB(2, ML)
      DIMENSION IWTPBF(3, MP), IWTLBF(3, ML), IWTSBF(3, ML)
      DIMENSION LINKPB(2, MP), LINKLB(2, ML), LINKSB(2, ML)

      DIMENSION LIST(NNUID), LISTN(NNUID), NUID(NNUID), XN(NPNODE)
      DIMENSION YN(NPNODE), NXK(NNXK, NPELEM), MAT(NPELEM)
      DIMENSION KXN(NNXK, MAXKXN), LA(NPNODE), LB(NPNODE)
      DIMENSION IHERE(NNUID), ILIST(MXLPS), XLIST(MXLPS)

      DIMENSION LSTNBC(NPNBC), LSTSBC(NPSBC), NSIDEN(NPSBC)
      DIMENSION NNFLG(MXNFLG), NNLEN(MXNFLG), NNPTR(MXNFLG)
      DIMENSION NSFLG(MXSFLG), NSLEN(MXSFLG), NSPTR(MXSFLG)
      DIMENSION WTHOLD(NPWTS), WTNODE(NPNBC), WTSIDE(NPSBC)
      DIMENSION NVLEN(MXSFLG), NVPTR(MXSFLG), MATMAP(3, NPREGN)

      DIMENSION KLIST(20)

      LOGICAL OPTIM, ERR, NOROOM, ALL, THREE, EIGHT, NINE
      LOGICAL ITSOK, ISBARS

      CHARACTER*80 NUMBER(MSC)

C  HEADER

      CALL MESAGE (' ')
      CALL MESAGE ('NUMBERING OF GENERATED OUTPUT BEGUN')
      IF (OPTIM) CALL MESAGE ('  -- OPTIMIZATION IS ENABLED --')

C  READ THE MESH TAPE

      CALL RDMESH (NPNODE, NPELEM, NPNBC, NPSBC, NPREGN, MS, MR, NNUID,
     &   NNXK, IUNIT, NNN, KKK, IPART, LSTNBC, LSTSBC, NUID, XN, YN,
     &   NXK, MAT, MATMAP, NUMMAT, ISIDE, NLPS, IFLINE, ILLIST, LINKS,
     &   LINKR, IMAT, LINKB, JMAT, NNNBC, NNSBC, ERR)
      IF (ERR) THEN
         CALL MESAGE ('** NUMBERING ABORT **')
         RETURN
      END IF

C  SORT NODE LIST INTO INCREASING NUID-S

      DO 100 I = 1, NNN
         LISTN(I) = NUID(I)
  100 CONTINUE
      IF (OPTIM) THEN
         CALL NODORD (NPNODE, XN, YN, LISTN, NUID, NNN)
      ELSE
         DO 110 I = 1, NNN
            LIST(I) = I
  110    CONTINUE
         CALL SORT (NNN, NUID, LIST)
      END IF

C  CONVERT REFERENCES TO NUID-S TO REFERENCES TO
C  SEQUENCE NUMBERS

      DO 130 I = 1, 4
         DO 120 K = 1, KKK
            IF (NXK(I, K) .GT. 0) THEN
               NEW = INDX(NNN, NUID, NXK(I, K))
               IF (NEW .EQ. 0) THEN
                  CALL MESAGE ('ERROR SORTING NUMBERING DATA -- NODE')
                  CALL MESAGE ('*** NO MESH SAVED ***')
                  KKK = 0
                  RETURN
               END IF
               IF (OPTIM) THEN
                  NXK(I, K) = NEW
               ELSE
                  NXK(I, K) = LIST(NEW)
               END IF
            END IF
  120    CONTINUE
  130 CONTINUE

      IF (NNNBC .GT. 0) THEN
         DO 140 I = 1, NNNBC
            IF (LSTNBC(I) .GT. 0) THEN
               NEW = INDX(NNN, NUID, LSTNBC(I))
               IF (NEW .EQ. 0) THEN
                  CALL MESAGE ('ERROR SORTING NUMBERING DATA -- NBC')
                  CALL MESAGE ('*** NO MESH SAVED ***')
                  KKK = 0
                  RETURN
               END IF
               IF (OPTIM) THEN
                  LSTNBC(I) = NEW
               ELSE
                  LSTNBC(I) = LIST(NEW)
               END IF
            END IF
  140    CONTINUE
      END IF

C  BUILD KXN ARRAY

      NUMKXN = NNN
      DO 160 I = 1, NNXK
         DO 150 J = 1, MAXKXN
            KXN(I, J) = 0
  150    CONTINUE
  160 CONTINUE

      DO 180 I = 1, 4
         DO 170 K = 1, KKK
            IF (NXK(I, K) .GT. 0) THEN
               CALL KXNADD (MAXKXN, NNXK, KXN, NUMKXN, K, NXK(I, K),
     &            ERR)
               IF (ERR) THEN
                  CALL MESAGE ('** NUMBERING ABORT **')
                  KKK = 0
                  RETURN
               END IF
            END IF
  170    CONTINUE
  180 CONTINUE
      IF (OPTIM) THEN

C  GET STARTING LIST FOR CUTHILL-MCKEE PROCESS

         IF (KCRD .GT. 0) THEN
            CALL GNLIST (NPNODE, NNUID, MSC, NPNODE, NPELEM, MAXKXN,
     &         NNXK, KXN, NXK, NUID, XN, YN, LIST, NUML, NUMBER, KCRD,
     &         NNN, ERR, NOROOM)
            IF (NOROOM) THEN
               CALL MESAGE ('TOO MANY NODES IN STARTING LIST')
               CALL MESAGE ('*** NO MESH SAVED ***')
               CALL MESAGE ('** NUMBERING ABORT **')
               KKK = 0
               RETURN
            ELSE IF (ERR) THEN
               CALL MESAGE ('ERROR GENERATING STARTING LIST')
               CALL MESAGE ('*** NO MESH SAVED ***')
               CALL MESAGE ('** NUMBERING ABORT **')
               KKK = 0
               RETURN
            END IF
         ELSE
            CALL MESAGE ('NO CARDS AVAILABLE FOR STARTING LIST')
            CALL MESAGE ('FIRST NODE IN THE OUTPUT USED')
            NUML = 1
            LIST(1) = 1
         END IF

C  INITIALIZE LISTS

         DO 190 I = 1, NNN
            LISTN(I) = I
  190    CONTINUE

C  USE LISTN AS A CHECK ON WHETHER THE NODE HAS BEEN USED (NEGATED)

         DO 200 I = 1, NUML
            LA(I) = LIST(I)
            NODE = LIST(I)
            LISTN(NODE) = -LISTN(NODE)
  200    CONTINUE

         NUMA = NUML

C  CREATE LIST OF NEW NODES CONNECTED TO LIST A

  210    CONTINUE
         NUMB = 0
         DO 230 N = 1, NUMA
            ALL = .TRUE.
            CALL GETNXN (NPNODE, NPELEM, MAXKXN, NNXK, KXN, NXK, LISTN,
     &         LA(N), NLIST, NUMB1, ALL, ERR)
            IF (ERR) THEN
               CALL MESAGE ('** NUMBERING ABORT **')
               KKK = 0
               RETURN
            END IF
            IF (NUMB1 .GT. 0) THEN
               IF ((NUMB + NUMB1) .GT. NPNODE) THEN
                  CALL MESAGE ('LIST B HAS OVERFLOWED')
                  CALL MESAGE ('*** NO MESH SAVED  ***')
                  CALL MESAGE ('** NUMBERING ABORT **')
                  KKK = 0
                  RETURN
               END IF
               DO 220 I = 1, NUMB1
                  NODE = NLIST(I)
                  NUMB = NUMB + 1
                  LB(NUMB) = NODE
                  LISTN(NODE) = -LISTN(NODE)
  220          CONTINUE
            END IF
  230    CONTINUE
         IF (NUMB .GT. 0) THEN

C  INCLUDE LIST B INTO FULL LIST
C  ALSO TRANSFER LIST B TO LIST A

            DO 240 I = 1, NUMB
               NUML = NUML + 1
               LIST(NUML) = LB(I)
               LA(I) = LB(I)
  240       CONTINUE
            NUMA = NUMB

C  CHECK FOR CONVERGENCE

            IF (NUML .LT. NNN) GO TO 210

C  PROCESS HAS CONVERGED
C  CHECK IF ALL NODES WERE COVERED

         ELSE IF (NUML .LT. NNN) THEN

            DO 250 I = 1, NNN
               IF (LISTN(I) .GT. 0) THEN

C  START THE LIST AGAIN WITH THE MISSED NODE

                  CALL MESAGE ('A DISCONTINUITY (SLIDE LINE) IN THE '//
     &               'BODY HAS BEEN FOUND')
                  CALL MESAGE
     &               ('SEPARATE PART NUMBERING WILL ALSO BE OPTIMIZED')
                  NUMA = 1
                  LA(1) = I
                  GO TO 210
               END IF
  250       CONTINUE

C  DEFINITE ERROR IN THE NUMBERING PROCESS

            CALL MESAGE ('ALL NODES COULD NOT BE FOUND TO NUMBER')
            CALL MESAGE ('       *** NO MESH SAVED ***')
            CALL MESAGE ('       ** NUMBERING ABORT **')
            KKK = 0
            RETURN
         END IF

C  PREPARE TO PUT NODE LIST INTO NEWLY DETERMINED ORDER
C  LISTN BECOMES THE POINTER FROM THE OLD NUMBER TO THE NEW

         DO 260 I = 1, NNN
            J = LIST(I)
            LISTN(J) = I
  260    CONTINUE

C  CONVERT NODE NUMBERS TO NEW NODE ORDER BY REDOING THE NXK ARRAY

         DO 280 I = 1, 4
            DO 270 K = 1, KKK
               J = NXK(I, K)
               IF (J .GT. 0) NXK(I, K) = LISTN(J)
  270       CONTINUE
  280    CONTINUE

         IF (NNNBC .GT. 0) THEN
            DO 290 I = 1, NNNBC
               IF (LSTNBC(I) .GT. 0) THEN
                  J = LSTNBC(I)
                  LSTNBC(I) = LISTN(J)
               END IF
  290       CONTINUE
         END IF

C  PUT NODE LIST INTO NEW ORDER

         CALL NODORD (NPNODE, XN, YN, LISTN, NUID, NNN)

C  REBUILD KXN ARRAY

         NUMKXN = NNN
         DO 310 I = 1, NNXK
            DO 300 J = 1, MAXKXN
               KXN(I, J) = 0
  300       CONTINUE
  310    CONTINUE

         DO 330 I = 1, 4
            DO 320 K = 1, KKK
               IF (NXK(I, K) .GT. 0) THEN
                  CALL KXNADD (MAXKXN, NNXK, KXN, NUMKXN, K, NXK(I, K),
     &               ERR)
                  IF (ERR) THEN
                     CALL MESAGE ('** NUMBERING ABORT **')
                     KKK = 0
                     RETURN
                  END IF
               END IF
  320       CONTINUE
  330    CONTINUE

C  PUT ELEMENT NUMBERING INTO NEW ORDER USING LA AS TEMPORARY STORAGE

         DO 340 I = 1, KKK
            LA(I) = 0
            LB(I) = 0
  340    CONTINUE
         KOUNT = 0
         DO 360 I = 1, NNN
            CALL GETKXN (NPNODE, MAXKXN, NNXK, KXN, NUID, I, KLIST,
     &         NUMK, ERR)
            DO 350 J = 1, NUMK
               IF (LB(KLIST(J)) .EQ. 0) THEN
                  KOUNT = KOUNT + 1
                  LA(KOUNT) = KLIST(J)
                  LB(KLIST(J)) = 1
               END IF
  350       CONTINUE
  360    CONTINUE

C  END OF OPTIMIZATION

      ELSE
         DO 370 I = 1, NNN
            NUID(I) = LISTN(I)
  370    CONTINUE
      END IF

C  STICK LSTNBC INTO LISTN AS A WORK ARRAY FOR SORTING NODAL BOUNDARY
C  CONDITIONS LISTS

      IF (NNNBC .GT. 0) THEN
         DO 380 I = 1, NNNBC
            LISTN(I) = LSTNBC(I)
  380    CONTINUE

C  SORT THROUGH LSTNBC AND RECREATE IT IN PLACE
C  USING LISTN AS THE ARRAY TO TAKE LSTNBC OVER
C  AND IHERE AS A WORK ARRAY
C  (LSTNBC NOW BECOMES THE NODES ARRAY FOR THE
C  GENESIS DATA BASE)

         CALL SRTNBC (MXNFLG, NPNBC, NNN, NNFLG, NNLEN, NNPTR, LSTNBC,
     &      LISTN, IHERE, NNNBC, NBCNOD, NNLIST)
      ELSE
         NNLIST = 0
         NBCNOD = 0
      END IF

C  SORT THROUGH LSTSBC AND RECREATE IT IN PLACE
C  USING LISTN AS THE ARRAY TO TAKE LSTSBC OVER
C  AND KXN AS A WORK ARRAY
C  (LSTSBC NOW BECOMES THE NELEMS ARRAY FOR THE
C  GENESIS DATA BASE)

      IF (NNSBC .GT. 0) THEN
         DO 390 I = 1, NNSBC
            LISTN(I) = LSTSBC(I)
  390    CONTINUE

         CALL SRTSBC (MXSFLG, NPSBC, NPELEM, NNXK, NXK, NSFLG, NSLEN,
     &      NSPTR, NVLEN, NVPTR, LISTN, LSTSBC, NSIDEN, IHERE, NNSBC,
     &      NSLIST, NVLIST, NBCSID)

      ELSE
         NBCSID = 0
         NSLIST = 0
         NVLIST = 0
      END IF

C  PUT WEIGHTS ON FLAGGED NODES AS NEEDED
C  USE THE IHERE ARRAY AS A WORK ARRAY

      CALL ADDWT (NNUID, NNXK, MAXKXN, NPNODE, NPELEM, MXLPS, MP, ML,
     &   MS, NPNBC, NPSBC, MXNFLG, MXSFLG, NPWTS, COOR, ILINE, LTYPE,
     &   LCON, ISIDE, NLPS, IFLINE, ILLIST, LINKP, LINKL, LINKS, IPBF,
     &   NPPF, IFPB, LISTPB, IWTPBF, ILBF, NLPF, IFLB, LISTLB, IWTLBF,
     &   ISBF, NSPF, IFSB, LISTSB, IWTSBF, LINKPB, LINKLB, LINKSB, XN,
     &   YN, NUID, NXK, KXN, LSTNBC, NNFLG, NNPTR, NNLEN, NSFLG, NVPTR,
     &   NVLEN, NSIDEN, WTNODE, WTSIDE, WTHOLD, NBCNOD, NNLIST, NBCSID,
     &   NSLIST, NVLIST, ILIST, XLIST)

C  SORT NUMBERS ACCORDING TO MATERIAL TYPE
C  USE KXN AS A WORK ARRAY

      DO 410 J = 1, 4
         DO 400 I = 1, KKK
            KXN(J, I) = NXK(J, I)
  400    CONTINUE
  410 CONTINUE

C  SET UP THE MATERIAL MAPPING ARRAY
C      MATMAP(1, I) = THE MATERIAL ID FOR THE I'TH BLOCK
C      MATMAP(2, I) = THE FIRST ELEMENT IN THE I'TH BLOCK
C      MATMAP(3, I) = THE LAST ELEMENT IN THE I'TH BLOCK

      KOUNT = 1
      DO 440 I = 1, NUMMAT
         KMAT = MATMAP(1, I)
         MATMAP(2, I) = KOUNT
         DO 430 J = 1, KKK
            IF (MAT(J) .EQ. KMAT) THEN
               LISTN(J) = KOUNT
               LIST(KOUNT) = J
               DO 420 K = 1, 4
                  NXK(K, KOUNT) = KXN(K, J)
  420          CONTINUE
               KOUNT = KOUNT + 1
            END IF
  430    CONTINUE
         MATMAP(3, I) = KOUNT - 1
  440 CONTINUE
      IF (KOUNT - 1 .NE. KKK) THEN
         CALL MESAGE ('ALL ELEMENTS DID NOT HAVE AN ELEMENT ID')
         CALL MESAGE ('MESH NUMBERING ABORTED')
         KKK = 0
         RETURN
      END IF

C  REDO THE REGION POINTER ARRAY

      DO 450 I = 1, NPREGN
         IPART(2, I) = LISTN(IPART(2, I))
         IPART(3, I) = LISTN(IPART(3, I))
  450 CONTINUE

C  REDO THE MATERIAL ARRAY

      DO 470 I = 1, NUMMAT
         DO 460 J = MATMAP(2, I), MATMAP(3, I)
            MAT(J) = MATMAP(1, I)
  460    CONTINUE
  470 CONTINUE

C  REDO THE MAPPING ARRAY IF OPTIMIZING RENUMBERING HAS BEEN DONE

      IF (OPTIM) THEN
         DO 480 I = 1, KKK
            LA(I) = LISTN(LA(I))
  480    CONTINUE
      END IF

C  REDO THE ELEMENT SIDE BOUNDARY LISTING WITH THE CURRENT ELEMENT NO.

      DO 490 I = 1, NSLIST
         LSTSBC(I) = LISTN(LSTSBC(I))
  490 CONTINUE

C  STORE THE LISTN POINTER SYSTEM FOR NOW

      DO 500 I = 1, KKK
         LIST(I) = LISTN(I)
  500 CONTINUE

C  ADD THE MID-SIDE NODES IF EIGHT OR NINE NODE QUADS ARE WANTED
C  OR IF THREE NODE BARS ARE WANTED

      IF ((EIGHT) .OR. (NINE) .OR. (THREE)) THEN

C  FLAG ALL ELEMENT SIDES ONLY ONCE (NO SHARED SIDE FLAGGED)

         CALL NXKBDY (NNXK * MAXKXN, NNXK, NPELEM, NXK, KKK, KXN,
     &      THREE, EIGHT, NINE)

C  CREATE THE MIDSIDE NODES

         CALL MIDNOD (NPNODE, NNUID, NPELEM, NNXK, MP, ML, KKK, NNN,
     &      NALL, NL, NXK, NUID, XN, YN, LISTN, COOR, ILINE, LTYPE,
     &      LCON, LINKP, LINKL, THREE, EIGHT, NINE)

C  MODIFY THE IDENTIFIERS OF THE OLD NODES

         DO 510 I = 1, NNN
            LISTN(I) = I * 100000
  510    CONTINUE
         NNN = NALL

C  ORDER THE EXPANDED NODE LIST

         CALL NODORD (NPNODE, XN, YN, LISTN, NUID, NNN)

C  EXPAND THE CONNECTIVITY ARRAY TO INCLUDE THE MIDSIDE NODES
C  WHILE REPOSITIONING THE CORNER NODES INTO PROPER SEQUENCE

         DO 540 I = 1, KKK
            IF ((THREE) .AND. (NXK (3, I) .EQ. 0)) THEN
               ITSOK = .TRUE.
            ELSEIF ((NXK (3, I) .NE. 0) .AND. ((EIGHT) .OR. (NINE)))
     &         THEN
               ITSOK = .TRUE.
            ELSE
               ITSOK = .FALSE.
            ENDIF
            IF (ITSOK) THEN
               DO 520 J = 4, 1, -1
                  JJ = J + 1
                  IF (JJ .EQ. 5) JJ = 1
                  NODEA = IABS(NXK(J, I))
                  NODEB = IABS(NXK(JJ, I))

C  CHECK FOR 3 NODE BAR ELEMENTS

                  IF ((NODEA .GT. 0) .AND. (NODEB .GT. 0)) THEN
                     NXK (J * 2 - 1, I) = NODEA * 100000
                     NLO = MIN0(NODEA, NODEB)
                     NHI = MAX0(NODEA, NODEB)
                     NXK (J * 2, I) = NLO * 100000 + NHI
                  ELSE IF (NODEA .GT. 0) THEN
                     NXK (J * 2 - 1, I) = NODEA * 100000
                  END IF
  520          CONTINUE
            ELSE
               DO 530 J = 1, 4
                  NXK (J, I) = IABS (NXK (J,I)) * 100000
  530          CONTINUE
            ENDIF
  540    CONTINUE

C  GET THE LIST OF NODAL BOUNDARY FLAGS EXTENDED AND THE NEW IDENTIFIERS
C  IN PLACE (USE NUID AND WTHOLD AS WORK ARRAYS)

         IF (NBCNOD .GT. 0) THEN
            KOUNT = 0
            DO 570 I = 1, NBCNOD
               JEND = NNLEN(I) + NNPTR(I) - 1
               DO 560 J = NNPTR(I), JEND
                  KOUNT = KOUNT + 1
                  NUID(KOUNT) = LSTNBC(J) * 100000
                  WTHOLD(KOUNT) = WTNODE(J)
                  DO 550 K = J + 1, JEND
                     ILO = MIN0(LSTNBC(J), LSTNBC(K))
                     IHI = MAX0(LSTNBC(J), LSTNBC(K))
                     ITRY = (ILO * 100000) + IHI
                     IF (INDX(NNN, LISTN, ITRY) .GT. 0) THEN
                        KOUNT = KOUNT + 1
                        NNLEN(I) = NNLEN(I) + 1
                        NUID(KOUNT) = ITRY
                        WTHOLD(KOUNT) = (WTNODE(J) + WTNODE(K)) * .5
                     END IF
  550             CONTINUE
  560          CONTINUE
               IF (I .NE. 1) NNPTR(I) = NNPTR(I - 1) + NNLEN(I - 1)
  570       CONTINUE
            NNLIST = KOUNT
            DO 580 I = 1, NNLIST
               LSTNBC(I) = NUID(I)
               WTNODE(I) = WTHOLD(I)
  580       CONTINUE
         END IF

C  GET THE LIST OF SIDE BOUNDARY FLAGS EXTENDED AND THE NEW IDENTIFIERS
C  IN PLACE (USE NUID, WTHOLD, AND KXN AS WORK ARRAYS)

         IF (NBCSID .GT. 0) THEN
            KOUNT = 0
            KOUNT2 = 0
            DO 600 I = 1, NBCSID
               ISTART = NSPTR (I)
               IEND = NSPTR (I) + NSLEN (I) - 1
               J = NVPTR (I)
               JEND = NVPTR (I) + NVLEN (I) - 1
               DO 590 II = ISTART, IEND
                  KELEM = LSTSBC (II)
                  KOUNT = KOUNT + 1
                  KXN (1, KOUNT) = KELEM
                  K = J + 1
                  KOUNT2 = KOUNT2 + 1
                  NUID(KOUNT2) = NSIDEN(J) * 100000
                  WTHOLD(KOUNT2) = WTSIDE(J)

C  DO THE ADJUSTMENTS IF THE ELEMENT IS ONE THAT HAS BEEN EXPANDED

                  IF ( ((THREE) .AND. (NXK (4, KELEM) .EQ. 0)) .OR.
     &               ( ((EIGHT) .OR. (NINE)) .AND.
     &               (NXK (4, KELEM) .NE. 0) ) ) THEN
                     KOUNT = KOUNT + 1
                     KXN (1, KOUNT) = LSTSBC(II)
                     NSLEN (I) = NSLEN (I) + 1
                     ILO = MIN0(NSIDEN(J), NSIDEN(K))
                     IHI = MAX0(NSIDEN(J), NSIDEN(K))
                     ITRY = (ILO * 100000) + IHI
                     KOUNT2 = KOUNT2 + 1
                     NUID(KOUNT2) = ITRY
                     WTHOLD(KOUNT2) = (WTSIDE(J) + WTSIDE(K)) * .5
                     KOUNT2 = KOUNT2 + 1
                     NUID(KOUNT2) = ITRY
                     WTHOLD(KOUNT2) = (WTSIDE(J) + WTSIDE(K)) * .5
                     NVLEN(I) = NVLEN(I) + 2
                  ENDIF
                  KOUNT2 = KOUNT2 + 1
                  NUID(KOUNT2) = NSIDEN(K) * 100000
                  WTHOLD(KOUNT2) = WTSIDE(K)
                  J = J + 2
  590          CONTINUE
               IF (I .GT. 1) THEN
                  NSPTR(I) = NSPTR(I - 1) + NSLEN(I - 1)
                  NVPTR(I) = NVPTR(I - 1) + NVLEN(I - 1)
               ENDIF
  600       CONTINUE

C  TRANSFER THE ELEMENT BOUNDARIES BACK FROM THE WORK ARRAYS

            NSLIST = KOUNT
            DO 610 I = 1, NSLIST
               LSTSBC(I) = KXN(1, I)
  610       CONTINUE
            NVLIST = KOUNT2
            DO 620 I = 1, NVLIST
               NSIDEN(I) = NUID(I)
               WTSIDE(I) = WTHOLD(I)
  620       CONTINUE
         END IF

C  ADD A CENTER NODE TO THE NODE LIST IF NEEDED

         IF (NINE) THEN
            NOLD = NNN
            DO 630 I = 1, KKK

C  WATCH OUT FOR 3 NODE BAR ELEMENTS

               IF (.NOT. ISBARS .OR. NXK(4, I) .GT. 0) THEN
                  N2 = INDX(NOLD, LISTN, NXK(2, I))
                  N4 = INDX(NOLD, LISTN, NXK(4, I))
                  N6 = INDX(NOLD, LISTN, NXK(6, I))
                  N8 = INDX(NOLD, LISTN, NXK(8, I))
                  IF ((N2 .EQ. 0) .OR. (N4 .EQ. 0) .OR.
     &               (N6 .EQ. 0) .OR. (N8 .EQ. 0)) THEN
                     CALL MESAGE ('BAD LINK IN RENUM AT 8 NODE LIST')
                     CALL MESAGE ('NO MESH SAVED')
                     KKK = 0
                     RETURN
                  END IF
                  DIST1 = SQRT((XN(N2) - XN(N6)) **2
     &               + (YN(N2) - YN(N6)) **2)
                  DIST2 = SQRT((XN(N8) - XN(N4)) **2
     &               + (YN(N8) - YN(N4)) **2)
                  NNN = NNN + 1
                  IF (DIST1 .LT. DIST2) THEN
                     XN(NNN) = .5 * (XN(N2) + XN(N6))
                     YN(NNN) = .5 * (YN(N2) + YN(N6))
                  ELSE
                     XN(NNN) = .5 * (XN(N4) + XN(N8))
                     YN(NNN) = .5 * (YN(N4) + YN(N8))
                  END IF
                  LISTN(NNN) = INT((MIN0(NXK(1, I), NXK(5, I))) +
     &               (MAX0(NXK(1, I), NXK(5, I)) * .0001))
                  NXK(9, I) = LISTN(NNN)
                ELSE
                  write (*,*) 'Element ', i, ' is a 3-node bar?'
               END IF
  630       CONTINUE

C  NOW, ORDER THE EXPANDED NODE LIST AGAIN

            CALL NODORD (NPNODE, XN, YN, LISTN, NUID, NNN)
            IEND = 9
         ELSEIF (EIGHT) THEN
            IEND = 8
         ELSE
            IEND = 4
         END IF

C  NOW REPLACE THE NODE REFERENCES WITH AN EXPANDED ORDER NUMBER

C  FIRST FIX THE CONNECTIVITY (NXK ARRAY)

         DO 650 I = 1, IEND
            DO 640 K = 1, KKK

C  AGAIN, WATCH OUT FOR 3 NODE BAR ELEMENTS

               IF (.NOT. ISBARS .OR. NXK(I, K) .GT. 0) THEN
                  NEW = INDX(NNN, LISTN, NXK(I, K))
                  IF (NEW .EQ. 0) THEN
                     CALL MESAGE ('BAD LINK IN RENUM 8 NODE LIST')
                     CALL MESAGE ('NO MESH SAVED')
                     KKK = 0
                     RETURN
                  END IF
                  NXK(I, K) = NEW
               END IF
  640       CONTINUE
  650    CONTINUE

C  NOW FIX THE NODE BOUNDARY FLAGS

         IF (NBCNOD .GT. 0) THEN
            DO 660 I = 1, NNLIST
               NEW = INDX(NNN, LISTN, LSTNBC(I))
               IF (NEW .EQ. 0) THEN
                  CALL MESAGE ('BAD LINK IN RENUM AT 8 NODE NBC')
                  CALL MESAGE ('NO MESH SAVED')
                  KKK = 0
                  RETURN
               END IF
               LSTNBC(I) = NEW
  660       CONTINUE
         END IF

C  NOW FIX THE SIDE BOUNDARY FLAGS

         IF (NBCSID .GT. 0) THEN
            DO 670 I = 1, NVLIST
               NEW = INDX(NNN, LISTN, NSIDEN(I))
               IF (NEW .EQ. 0) THEN
                  CALL MESAGE ('BAD LINK IN RENUM AT 8 NODE SBC')
                  CALL MESAGE ('NO MESH SAVED')
                  KKK = 0
                  RETURN
               END IF
               NSIDEN(I) = NEW
  670       CONTINUE
         END IF
      END IF

C  RENUMBERING COMPLETED

      CALL MESAGE (' ')
      CALL MESAGE ('**************************************************')
      CALL MESAGE ('**           MESH PROCESSING COMPLETED          **')
      IF (ISBARS .AND. NINE) THEN
         CALL MESAGE ('**            THREE NODE BARS OUTPUT         '//
     &      '   **')
      ENDIF
      IF (NINE) THEN
         CALL MESAGE ('**            NINE NODE QUADS OUTPUT         '//
     &      '   **')
      ELSE IF (EIGHT) THEN
         CALL MESAGE ('**           EIGHT NODE QUADS OUTPUT         '//
     &      '   **')
      END IF
      IF (OPTIM) THEN
         CALL MESAGE
     &      ('**   WITH NODE AND ELEMENT NUMBERING OPTIMIZED  **')

C  FIND LARGEST NODE DIFFERENCE FOR AN ELEMENT

         LWID = 0
         DO 680 K = 1, KKK
            N1 = NXK(1, K)
            N2 = NXK(2, K)
            N3 = NXK(3, K)
            N4 = NXK(4, K)
            IF ((N4 .GT. 0) .AND. ((EIGHT) .OR. (NINE))) THEN
               N5 = NXK(5, K)
               N6 = NXK(6, K)
               N7 = NXK(7, K)
               N8 = NXK(8, K)
               IF (NINE) THEN
                  N9 = NXK(9, K)
                  NLO = MIN0(N1, N2, N3, N4, N5, N6, N7, N8, N9)
                  NHI = MAX0(N1, N2, N3, N4, N5, N6, N7, N8, N9)
               ELSE
                  NLO = MIN0(N1, N2, N3, N4, N5, N6, N7, N8)
                  NHI = MAX0(N1, N2, N3, N4, N5, N6, N7, N8)
               END IF
            ELSE
               IF (N3 .LE. 0) N3 = N1
               IF (N4 .LE. 0) N4 = N2
               NLO = MIN0(N1, N2, N3, N4)
               NHI = MAX0(N1, N2, N3, N4)
            END IF
            LWID = MAX0(LWID, NHI - NLO)
  680    CONTINUE
         WRITE(*, 10000) LWID
      END IF
      WRITE(*, 10010) NNN, KKK, NUMMAT
      CALL MESAGE ('**************************************************')

C  RESTORE THE NUID ARRAY AS A POINTER ARRAY OF OLD TO NEW ELEMENTS
C  (MAPDXG ARRAY)

      IF (OPTIM) THEN
         DO 690 I = 1, KKK
            NUID(I) = LA(I)
            LIST(NUID(I)) = I
  690    CONTINUE
      ELSE
         DO 700 I = 1, KKK
            NUID(I) = I
            LIST(I) = I
  700    CONTINUE
      END IF

      RETURN

10000 FORMAT(' **   LARGEST NODE DIFFERENCE PER ELEMENT:', I6, '  **')
10010 FORMAT(' ** NODES:', I6, '; ELEMENTS:', I6, '; MATERIALS:', I3,
     &   ' **')

      END
