C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ADDWT (NNUID, NNXK, MAXKXN, NPNODE, NPELEM, MXLPS, MP,
     &   ML, MS, NPNBC, NPSBC, MXNFLG, MXSFLG, NPWTS, COOR, ILINE,
     &   LTYPE, LCON, ISIDE, NLPS, IFLINE, ILLIST, LINKP, LINKL, LINKS,
     &   IPBF, NPPF, IFPB, LISTPB, IWTPBF, ILBF, NLPF, IFLB, LISTLB,
     &   IWTLBF, ISBF, NSPF, IFSB, LISTSB, IWTSBF, LINKPB, LINKLB,
     &   LINKSB, XN, YN, NUID, NXK, KXN, LSTNBC, NNFLG, NNPTR, NNLEN,
     &   NSFLG, NVPTR, NVLEN, NSIDEN, WTNODE, WTSIDE, WTHOLD, NBCNOD,
     &   NNLIST, NBCSID, NSLIST, NVLIST, ILIST, XLIST)
C***********************************************************************

C  SUBROUTINE ADDWT = ADDS THE WEIGHTING FACTORS TO ANY NODES WITH
C                     FLAGS CONTAINING WEIGHTS

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     RENUM = NUMBERS QMESH OUTPUT,  AND RENUMBERS AS NEEDED FOR
C             OPTIMIZATION

C***********************************************************************

      DIMENSION COOR (2, MP), ILINE (ML), LTYPE (ML), LCON (3, ML)
      DIMENSION ISIDE (MS), NLPS (MS), IFLINE (MS), ILLIST (MS*3)
      DIMENSION LINKP (2, MP), LINKL (2, ML), LINKS (2, MS)
      DIMENSION IPBF (MP), NPPF (MP), IFPB (MP), LISTPB (2, MP),
     &   IWTPBF (3, MP)
      DIMENSION ILBF (ML), NLPF (ML), IFLB (ML), LISTLB (2, ML),
     &   IWTLBF (3, ML)
      DIMENSION ISBF (ML), NSPF (ML), IFSB (ML), LISTSB (2, ML),
     &   IWTSBF (3, ML)
      DIMENSION LINKPB (2, MP), LINKLB (2, ML), LINKSB (2, ML)

      DIMENSION NUID (NNUID), NXK (NNXK, NPELEM), KXN (NNXK, MAXKXN)
      DIMENSION XN (NPNODE), YN (NPNODE), ILIST (MXLPS), XLIST (MXLPS)

      DIMENSION LSTNBC (NPNBC), NSIDEN (NPSBC), WTHOLD (NPWTS)
      DIMENSION NNFLG (MXNFLG), NNLEN (MXNFLG), NNPTR (MXNFLG),
     &   WTNODE (NPNBC)
      DIMENSION NSFLG (MXSFLG), NVLEN (MXSFLG), NVPTR (MXSFLG),
     &   WTSIDE (NPSBC)

      LOGICAL ADDLNK, ERR, ISPNT

      ADDLNK = .FALSE.
      IPNTR  = 0

C FIRST FLAG ALL WEIGHT ARRAYS TO -1.0 TO KNOW WHICH REMAIN DEFAULTED

      DO 100 I = 1, NNLIST
         WTNODE (I) = -1.0
         WTHOLD (I) = -1.0
  100 CONTINUE
      DO 110 I = 1, NVLIST
         WTSIDE (I) = -1.0
  110 CONTINUE

C  NOW CHECK ALL POINT FLAGS FOR WEIGHTS AND APPLY THE POINT
C  Y VALUE AS THE WEIGHT FOR THE NODE AT THE BEGINNING POINT

      ISPNT = .TRUE.
      DO 120 I = 1, NBCNOD
         CALL LTSORT (MP, LINKPB, NNFLG (I), IPNTR, ADDLNK)
         IF (IPNTR .GT. 0) THEN
            IF (IWTPBF (1, IPNTR) .GT. 0) THEN
               IPOINT = IWTPBF (1, IPNTR)
               JPOINT = IWTPBF (2, IPNTR)
               CALL CHKWT (MP, ML, MS, NNLIST, NBCNOD, NNUID, MXLPS,
     &            LINKP, LINKL, LINKS, NUID, NNFLG, NNLEN, NNPTR,
     &            LSTNBC, IPBF (IPNTR), LTYPE, LCON, NLPS, IFLINE,
     &            ILLIST, COOR, JPOINT, IPOINT, ILOC, JLOC, NIX, ILIST,
     &            XLIST, ADDLNK, ISPNT, ERR)
               IF (ERR) THEN
                  WTNODE (JLOC) = 1.0
               ELSE
                  CALL LTSORT (MP, LINKP, IPOINT, IPNTR, ADDLNK)
                  WTNODE (JLOC) = COOR (2, IPNTR)
                  WTHOLD (JLOC) = 1.
               ENDIF
            ENDIF
         ENDIF
  120 CONTINUE
      ISPNT = .FALSE.

C  NOW CHECK ALL LINE FLAGS FOR WEIGHTS AND APPLY THE APPROPRIATE
C  WEIGHT ALL ALONG CONTINUOUS NODES ON THE BOUNDARY.

      DO 160 I = 1, NBCNOD
         CALL LTSORT (ML, LINKLB, NNFLG (I), IPNTR, ADDLNK)
         IF (IPNTR .GT. 0) THEN
            IF (IWTLBF (1, IPNTR) .NE. 0) THEN
               WRITE (*, 10000)NNFLG (I)
               JSIDE = IWTLBF (1, IPNTR)
               JPOINT = IWTLBF (2, IPNTR)
               CALL CHKWT (MP, ML, MS, NNLIST, NBCNOD, NNUID, MXLPS,
     &            LINKP, LINKL, LINKS, NUID, NNFLG, NNLEN, NNPTR,
     &            LSTNBC, ILBF (IPNTR), LTYPE, LCON, NLPS, IFLINE,
     &            ILLIST, COOR, JPOINT, JSIDE, ILOC, JLOC, NIX, ILIST,
     &            XLIST, ADDLNK, ISPNT, ERR)
               IF (.NOT.ERR) THEN

C  LOOP UNTIL ALL THE NODES HAVE BEEN FOUND,
C  FIRST PUTTING THE ACCUMULATED LENGTH IN THE WTNODE ARRAY

                  NEWNOD = 0
                  ACCLEN = 0.
                  NODE = LSTNBC (JLOC)
                  IF (IWTLBF (3, IPNTR) .NE. 0) THEN
                     LINE1 = IWTLBF (3, IPNTR)
                  ELSE
                     LINE1 = LISTLB (1, IFLB (IPNTR))
                  ENDIF
  130             CONTINUE
                  WTNODE (JLOC) = ACCLEN
                  CALL BFNODE (NNLEN (ILOC), NNXK, NPNODE, NPELEM,
     &               MAXKXN, NNUID, NODE, NEWNOD, LSTNBC (NNPTR (ILOC)),
     &               KXN, NXK, NUID, JLOC, LINE1, ERR)
                  IF (ERR) THEN
                     WRITE (*, 10010)NNFLG (IPNTR)
                     GOTO 150
                  ENDIF
                  JLOC = JLOC+NNPTR (ILOC)-1
                  IF (NEWNOD .GT. 0) THEN
                     ACCLEN = ACCLEN+SQRT (( (XN (NODE)-XN (NEWNOD))
     &                  ** 2) + ((YN (NODE)-YN (NEWNOD)) ** 2))
                     NUID (NODE) = -ABS (NUID (NODE))
                     NODE = NEWNOD
                     NEWNOD = 0
                     GOTO 130
                  ENDIF

C  NOW CHANGE THE ACCUMULATED LENGTH TO A PERCENTAGE LENGTH
C  AND GET THE WEIGHTING FUNCTION

                  DO 140 J = NNPTR (ILOC), NNPTR (ILOC)+NNLEN (ILOC)-1
                     IF ((WTNODE (J) .GE. 0.).AND. (ACCLEN .NE. 0.))
     &                  THEN
                        WTHOLD (J) = 1.
                        WTNODE (J) = WTNODE (J)/ACCLEN
                        CALL GETWT (MP, ML, MXLPS, NIX, ILIST, XLIST,
     &                     ILINE, LCON, LTYPE, COOR, LINKP, WTNODE (J),
     &                     ADDLNK, ERR)
                        IF (ERR)GOTO 150
                     ENDIF
  140             CONTINUE
               ENDIF
            ENDIF
         ENDIF
  150    CONTINUE
  160 CONTINUE

C  NOW RESET NUIDS AND PUT ALL DEFAULTS TO 1.0

      DO 170 I = 1, NNLIST
         NUID (LSTNBC (I)) = ABS (NUID (LSTNBC (I)))
         IF (WTHOLD (I) .LT. 0.)WTNODE (I) = 1.0
  170 CONTINUE

C  NOW CHECK ALL SIDE FLAGS FOR WEIGHTS AND APPLY THE APPROPRIATE
C  WEIGHT ALL ALONG CONTINUOUS NODES ON THE BOUNDARY.

      DO 260 I = 1, NBCSID
         CALL LTSORT (MP, LINKSB, NSFLG (I), IPNTR, ADDLNK)
         IF (IPNTR .GT. 0) THEN
            IF (IWTSBF (1, IPNTR) .NE. 0) THEN
               WRITE (*, 10020)NSFLG (I)
               JSIDE = IWTSBF (1, IPNTR)
               JPOINT = IWTSBF (2, IPNTR)
               CALL CHKWT (MP, ML, MS, NVLIST, NBCSID, NNUID, MXLPS,
     &            LINKP, LINKL, LINKS, NUID, NSFLG, NVLEN, NVPTR,
     &            NSIDEN, ISBF (IPNTR), LTYPE, LCON, NLPS, IFLINE,
     &            ILLIST, COOR, JPOINT, JSIDE, ILOC, JLOC, NIX, ILIST,
     &            XLIST, ADDLNK, ISPNT, ERR)
               IF (ERR) THEN
                  DO 180 J = NVPTR (IPNTR), NVPTR (IPNTR) +
     &               NVLEN (IPNTR) + 1
                     WTSIDE (J) = 1.0
  180             CONTINUE
               ELSE

C  LOOP UNTIL ALL THE NODES HAVE BEEN FOUND,
C  FIRST PUTTING THE ACCUMULATED LENGTH IN THE WTSIDE ARRAY

                  NEWNOD = 0
                  J1 = NVPTR (ILOC)
                  J2 = NVPTR (ILOC)+NVLEN (ILOC)-1
                  ACCLEN = 0.
                  NODE = NSIDEN (JLOC)
                  IF (IWTSBF (3, IPNTR) .NE. 0) THEN
                     LINE1 = IWTSBF (3, IPNTR)
                  ELSE
                     LINE1 = LISTLB (1, IFSB (IPNTR))
                  ENDIF
  190             CONTINUE

C  PUT THIS ACCLEN FOR ALL OCCURRENCES OF NODE IN THE LIST

                  DO 200 J = J1, J2
                     IF (NSIDEN (J) .EQ. NODE)WTSIDE (J) = ACCLEN
  200             CONTINUE
                  CALL BFNODE (NVLEN (ILOC), NNXK, NPNODE, NPELEM,
     &               MAXKXN, NNUID, NODE, NEWNOD, NSIDEN (NVPTR (ILOC)),
     &               KXN, NXK, NUID, JLOC, LINE1, ERR)
                  IF (ERR) THEN
                     DO 210 J = NVPTR (IPNTR), NVPTR (IPNTR) +
     &                  NVLEN (IPNTR)+1
                        WTSIDE (J) = 1.0
  210                CONTINUE
                     WRITE (*, 10030)NSFLG (IPNTR)
                     GOTO 250
                  ENDIF
                  JLOC = JLOC+NNPTR (ILOC)-1
                  IF (NEWNOD .GT. 0) THEN
                     ACCLEN = ACCLEN + SQRT (( (XN (NODE)-XN (NEWNOD))
     &                  ** 2) + ((YN (NODE)-YN (NEWNOD)) ** 2))
                     NUID (NODE) = -ABS (NUID (NODE))
                     NODE = NEWNOD
                     NEWNOD = 0
                     GOTO 190
                  ENDIF

C  NOW CHANGE THE ACCUMULATED LENGTH TO A PERCENTAGE LENGTH
C  AND GET THE WEIGHTING FUNCTION

                  DO 220 J = NVPTR (ILOC), NVPTR (ILOC)+NVLEN (ILOC)-1
                     IF ((WTSIDE (J) .GE. 0.).AND. (ACCLEN .NE. 0.))
     &                  THEN
                        WTSIDE (J) = WTSIDE (J)/ACCLEN
                        CALL GETWT (MP, ML, MXLPS, NIX, ILIST, XLIST,
     &                     ILINE, LCON, LTYPE, COOR, LINKP, WTSIDE (J),
     &                     ADDLNK, ERR)
                        IF (ERR)GOTO 250
                     ELSE
                        WTSIDE (J) = 0.
                     ENDIF
  220             CONTINUE

C  NOW RESET NUIDS FROM THIS SIDE SET

                  DO 230 JI = 1, NVLIST
                     NUID (NSIDEN (JI)) = ABS (NUID (NSIDEN (JI)))
  230             CONTINUE

               ENDIF
            ELSE
               DO 240 J = NVPTR (I), NVPTR (I)+NVLEN (I)+1
                  WTSIDE (J) = 1.0
  240          CONTINUE
            ENDIF
         ENDIF
  250    CONTINUE
  260 CONTINUE

C  NOW RESET NUIDS

      DO 270 I = 1, NVLIST
         NUID (NSIDEN (I)) = ABS (NUID (NSIDEN (I)))
  270 CONTINUE

      RETURN

10000 FORMAT (/, ' WEIGHTING BEGUN FOR NODAL FLAG', I5)
10010 FORMAT (' NO WEIGHTING POSSIBLE FOR NODAL FLAG', I5)
10020 FORMAT (/, ' WEIGHTING BEGUN FOR ELEMENT FLAG', I5)
10030 FORMAT (' NO WEIGHTING POSSIBLE FOR ELEMENT FLAG', I5)

      END
