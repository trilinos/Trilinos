C $Id: rdgens.f,v 1.5 2001/02/26 16:56:51 gdsjaar Exp $
C $Log: rdgens.f,v $
C Revision 1.5  2001/02/26 16:56:51  gdsjaar
C Fixes to all 'mesh read' to work.  Problem was a write not inside an if(remesh) block.
C
C Upped version to 3.11
C
C Revision 1.4  2000/11/13 15:39:05  gdsjaar
C Cleaned up unused variables and labels.
C
C Removed some real to int conversion warnings.
C
C Revision 1.3  1998/07/14 18:19:52  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.2  1990/11/30 11:30:10  gdsjaar
C Rewrote indexing for reads and writes
C
c Revision 1.1.1.1  90/11/30  11:14:28  gdsjaar
c FASTQ Version 2.0X
c 
c Revision 1.1  90/11/30  11:14:27  gdsjaar
c Initial revision
c 
CC* FILE: [.MAIN]RDGENS.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE RDGENS (A, IA, NPNODE, NPELEM, MXNFLG, MXSFLG, NPNBC,
     &   NPSBC, NPREGN, IUNIT, NNN, KKK, NUMMAT,  MCOM, ICOM, JCOM,
     &   CIN, RIN, IIN, KIN, IDUMP, IPART, NODES, NELEMS, NNFLG,
     &   NNPTR, NNLEN, NSFLG, NSPTR, NSLEN, NVPTR, NVLEN, NSIDEN,
     &   MAPDXG, XN, YN, NXK, MAT, MAPGXD, MATMAP, WTNODE, WTSIDE,
     &   NAMEHI, NAMEGL, NAMENV, NAMEEV, MAXNAM, NBCNOD, NNLIST, NBCSID,
     &   NSLIST, NVLIST, TITLE, EIGHT, NINE, REMESH, ERR, END, AMESUR,
     &   LINKEG, LISTEG, BMESUR, CMESUR, DMESUR, ISEVOK, MLINK, NPNOLD,
     &   NPEOLD, REXMIN, REXMAX, REYMIN, REYMAX, IDIVIS, SIZMIN, EMAX,
     &   EMIN)
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/31/90
CC* MODIFICATION: ADDED ARGUMENTS TO PASS MINIMUM ELEMENT SIZE (SIZMIN)
CC**              AND GETSIZ PARAMETERS OF EMIN AND EMAX
C
C************************************************************************
C
C  SUBROUTINE RDGENS = READS EXODUS DATABASE MESH OUTPUT
C
C***********************************************************************
C
      DIMENSION A(1), IA(1)
      DIMENSION XN (NPNODE), YN (NPNODE), NXK (9, NPELEM), MAT (NPELEM)
      DIMENSION NODES (NPNBC), NELEMS (NPSBC), NSIDEN (NPSBC)
      DIMENSION NNFLG (MXNFLG), NNLEN (MXNFLG)
      DIMENSION NNPTR (MXNFLG), WTNODE (NPNBC)
      DIMENSION NSFLG (MXSFLG), NSLEN (MXSFLG)
      DIMENSION NSPTR (MXSFLG), WTSIDE (NPSBC)
      DIMENSION NVLEN (MXSFLG), NVPTR (MXSFLG), IPART (3, NPREGN)
      DIMENSION MAPDXG (NPNODE), MAPGXD (NPNODE), MATMAP (3, NPREGN)
      DIMENSION KIN(MCOM), IIN(MCOM), RIN(MCOM)
C
      DIMENSION ISEVOK (MAXNAM, NUMMAT)
      DIMENSION AMESUR (NPEOLD), BMESUR (NPNOLD)
      DIMENSION CMESUR (NPEOLD), DMESUR (NPNOLD)
      DIMENSION LINKEG(2, MLINK), LISTEG(4 * NPEOLD)
C
      CHARACTER*72 CIN(MCOM), TITLE
      CHARACTER*8 DUMQA(4)
      CHARACTER*80 DUMINF
      CHARACTER*8 NAMEHI(MAXNAM), NAMEGL(MAXNAM)
      CHARACTER*8 NAMENV(MAXNAM), NAMEEV(MAXNAM)
C
      LOGICAL EIGHT, NINE, ERR, END, REMESH
      LOGICAL FOUND1, FOUND2, SECOND, TERROR, MINSIZ
C
      integer lcon(9)
      data lcon /1,3,5,7,2,4,6,8,9/

      EIGHT = .FALSE.
      ERR = .TRUE.
      NINE = .FALSE.
      END = .TRUE.
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/17/90
CC* MODIFICATION: ADDED AN EXODUS WRITE FOR THE REMESH EXPANSION FACTOR.
C**               IT WILL BE REMOVED ONCE THE ROUTINE IS PROVEN.
C
C  READ THE NODE BLOCK
C
      READ (IUNIT, ERR = 270) (XN (I), I = 1, NNN),
     &   (YN (I), I = 1, NNN)
      READ (IUNIT, ERR = 270) (MAPDXG (I), I = 1, KKK)
      IF (REMESH) THEN
         IIUNIT = 21
         WRITE (IIUNIT, ERR = 270) (XN (I), I = 1, NNN),
     &      (YN (I), I = 1, NNN)
         WRITE (IIUNIT, ERR = 270) (MAPDXG (I), I = 1, KKK)
      ENDIF
C
C  READ THE ELEMENT BLOCKS
C
      KOUNT = 0
      DO 120 I = 1, NUMMAT
         READ (IUNIT, ERR = 270) MATMAP (1, I), NBLOCK, NNXK, IDUM
         IF (REMESH)
     &      WRITE (IIUNIT, ERR = 270) MATMAP (1, I), NBLOCK, NNXK, IDUM
         IF (NNXK .EQ. 8) EIGHT = .TRUE.
         IF (NNXK .EQ. 9) NINE = .TRUE.
         MATMAP (2, I) = KOUNT + 1
         KOUNT = KOUNT + NBLOCK
         MATMAP (3, I) = KOUNT
         IPART (1, I) = I
         IPART (2, I) = MATMAP (2, I)
         IPART (3, I) = MATMAP (3, I)
         IF (NNXK .EQ. 8) THEN
            read (iunit, err = 270) ((nxk(lcon(ii), k), ii=1, nnxk),
     $           k = matmap(2,i), matmap(3,i))
c$$$            READ (IUNIT, ERR = 270) ((NXK (1, K), NXK (3, K),
c$$$     &         NXK (5, K), NXK (7, K), NXK (2, K), NXK (4, K),
c$$$     &         NXK (6, K), NXK (8, K)),
c$$$     &         K = MATMAP (2, I), MATMAP (3, I))
         ELSEIF (NNXK .EQ. 9) THEN
            read (iunit, err = 270) ((nxk(lcon(ii), k), ii=1, nnxk),
     $           k = matmap(2,i), matmap(3,i))
c$$$            READ (IUNIT, ERR = 270) ((NXK (1, K), NXK (3, K),
c$$$     &         NXK (5, K), NXK (7, K), NXK (2, K), NXK (4, K),
c$$$     &         NXK (6, K), NXK (8, K), NXK (9, K)),
c$$$     &         K = MATMAP (2, I), MATMAP (3, I))
         ELSE
            READ (IUNIT, ERR = 270) ((NXK (J, K), J = 1, NNXK),
     &         K = MATMAP (2, I), MATMAP (3, I))
            IF (REMESH)
     &         WRITE (IIUNIT, ERR = 270) ((NXK (J, K), J = 1, NNXK),
     &         K = MATMAP (2, I), MATMAP (3, I))
         ENDIF
         READ (IUNIT, ERR = 270)
         IF (REMESH)
     *     WRITE (IIUNIT, ERR = 270) ((1., J = 1, IDUM), K = 1, KKK)
         DO 110 II = NNXK + 1, 9
            DO 100 J = MATMAP (2, I), MATMAP (3, I)
               NXK (II, J) = 0
  100       CONTINUE
  110    CONTINUE
  120 CONTINUE
C
C  READ THE NODAL BOUNDARY FLAGS
C
      READ (IUNIT, ERR = 270) (NNFLG (I), I = 1, NBCNOD)
      READ (IUNIT, ERR = 270) (NNLEN (I), I = 1, NBCNOD)
      READ (IUNIT, ERR = 270) (NNPTR (I), I = 1, NBCNOD)
      READ (IUNIT, ERR = 270) (NODES (I), I = 1, NNLIST)
      READ (IUNIT, ERR = 270) (WTNODE (I), I = 1, NNLIST)
C
C  READ THE SIDE BOUNDARY FLAGS
C
      READ (IUNIT, ERR = 270) (NSFLG (I), I = 1, NBCSID)
      READ (IUNIT, ERR = 270) (NSLEN (I), I = 1, NBCSID)
      READ (IUNIT, ERR = 270) (NVLEN (I), I = 1, NBCSID)
      READ (IUNIT, ERR = 270) (NSPTR (I), I = 1, NBCSID)
      READ (IUNIT, ERR = 270) (NVPTR (I), I = 1, NBCSID)
      READ (IUNIT, ERR = 270) (NELEMS (I), I = 1, NSLIST)
      READ (IUNIT, ERR = 270) (NSIDEN (I), I = 1, NVLIST)
      READ (IUNIT, ERR = 270) (WTSIDE (I), I = 1, NVLIST)
C
C  WRITE THE BOUNDARY FLAGS
C
      IF (REMESH) THEN
         WRITE (IIUNIT, ERR = 270) (NNFLG (I), I = 1, NBCNOD)
         WRITE (IIUNIT, ERR = 270) (NNLEN (I), I = 1, NBCNOD)
         WRITE (IIUNIT, ERR = 270) (NNPTR (I), I = 1, NBCNOD)
         WRITE (IIUNIT, ERR = 270) (NODES (I), I = 1, NNLIST)
         WRITE (IIUNIT, ERR = 270) (WTNODE (I), I = 1, NNLIST)
         WRITE (IIUNIT, ERR = 270) (NSFLG (I), I = 1, NBCSID)
         WRITE (IIUNIT, ERR = 270) (NSLEN (I), I = 1, NBCSID)
         WRITE (IIUNIT, ERR = 270) (NVLEN (I), I = 1, NBCSID)
         WRITE (IIUNIT, ERR = 270) (NSPTR (I), I = 1, NBCSID)
         WRITE (IIUNIT, ERR = 270) (NVPTR (I), I = 1, NBCSID)
         WRITE (IIUNIT, ERR = 270) (NELEMS (I), I = 1, NSLIST)
         WRITE (IIUNIT, ERR = 270) (NSIDEN (I), I = 1, NVLIST)
         WRITE (IIUNIT, ERR = 270) (WTSIDE (I), I = 1, NVLIST)
      ENDIF
C
C  READ THE RESULTS DATA IF REMESHING A PROBLEM
C
      IF (REMESH) THEN
C
C  READ IN THE QA, INFORMATION LINES AND
C  COORDINATE AND ELEMENT TYPE NAMES
C
         READ (IUNIT, ERR = 270) NQAREC
         WRITE (IIUNIT, ERR = 270) NQAREC
         DO 130 I = 1, NQAREC
            READ (IUNIT, ERR = 270) (DUMQA (J), J = 1, 4)
            WRITE (IIUNIT, ERR = 270) (DUMQA (J), J = 1, 4)
  130    CONTINUE
         READ (IUNIT, ERR = 270) NINFO
         WRITE (IIUNIT, ERR = 270) NINFO
         DO 140 I = 1, NINFO
            READ (IUNIT, ERR = 270) DUMINF
            WRITE (IIUNIT, ERR = 270) DUMINF
  140    CONTINUE
         READ (IUNIT, ERR = 270) (DUMQA (J), J = 1, 2)
         READ (IUNIT, ERR = 270) (DUMQA (1), J = 1, NUMMAT)
         WRITE (IIUNIT, ERR = 270) (DUMQA (J), J = 1, 2)
         WRITE (IIUNIT, ERR = 270) (DUMQA (1), J = 1, NUMMAT)
C
C  READ IN THE RESULTS VARIABLE NAMES
C
         READ (IUNIT, ERR = 270) NVARHI, NVARGL, NVARNP, NVAREL
         WRITE (IIUNIT, ERR = 270) 0, 0, 2, 0
         IF (MAX0 (NVARHI, NVARGL, NVARNP, NVAREL) .GT. MAXNAM) THEN
            CALL MESAGE ('** NOT ENOUGH ROOM FOR THE VARIABLE NAMES **')
            GOTO 270
         ENDIF
         READ (IUNIT, ERR= 270)
     &      (NAMEHI (I), I = 1, NVARHI),
     &      (NAMEGL (I), I = 1, NVARGL),
     &      (NAMENV (I), I = 1, NVARNP),
     &      (NAMEEV (I), I = 1, NVAREL)
         NAMENV(1) = 'ERR     '
         NAMENV(2) = 'EXP     '
         WRITE (IIUNIT, ERR= 270)
     &      (NAMEHI (I), I = 1, 0),
     &      (NAMEGL (I), I = 1, 0),
     &      (NAMENV (I), I = 1, 2),
     &      (NAMEEV (I), I = 1, 0)
C
C  GET THE NODAL VARIABLE OR ELEMENT VARIABLE NAME TO BE PROCESSED
C
         SECOND = .FALSE.
         IF (ICOM.GT.JCOM) THEN
            CALL MESAGE (' ')
            CALL MESAGE ('HERE ARE THE AVAILABLE NODAL VARIABLES: ')
            DO 160 I = 1, NVARNP
               WRITE (*, 10000) I, NAMENV (I)
  160       CONTINUE
            CALL MESAGE ('HERE ARE THE AVAILABLE ELEMENT VARIABLES: ')
            DO 170 I = 1, NVAREL
               WRITE (*, 10000) I + NVARNP, NAMEEV (I)
  170       CONTINUE
            IF (SECOND) THEN
               CALL FREFLD (IZ, IZ, 'ENTER SECOND ERROR ESTIMATOR '//
     &            'VARAIABLE NUMBER: ', MCOM, IOSTAT, JCOM, KIN, CIN,
     &            IIN, RIN)
            ELSE
               CALL FREFLD (IZ, IZ, 'ENTER ERROR ESTIMATOR VARAIABLE'//
     &            '  NUMBER: ', MCOM, IOSTAT, JCOM, KIN, CIN, IIN, RIN)
            ENDIF
            ICOM = 1
         END IF
         IVAR = IIN (ICOM)
         ICOM = ICOM + 1
         IF ((IVAR .LT. 1) .OR. (IVAR .GT. NVARNP + NVAREL)) THEN
            CALL MESAGE ('** NODAL VARIABLE NUMBER CHOSEN IS NOT '//
     &         'VALID **')
            GOTO 270
         ENDIF
C
C  CHECK TO SEE IF A SECOND VARIABLE IS DEISRED
C
         IF (.NOT. SECOND) THEN
            IVAR1 = IVAR
            IVAR2 = 0
         ELSE
            IVAR2 = IVAR
         ENDIF
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/12/90
CC* MODIFICATION: ADDED TARGET ERROR (TERROR & TERR) INQUIRIES
CC**              FOR ERROR ESTIMATION ROUTINES INCLUDING A MODIFICATION
CC**              TO THE CALL TO MIXVAR.
C
C  CHECK TO SEE IF A TARGET ERROR IS NEEDED
C
         CALL INTRUP ('USE A TARGET ERROR', TERROR, MCOM, ICOM, JCOM,
     &      CIN, IIN, RIN, KIN)
         IF (TERROR) THEN
            IF (ICOM .GT. JCOM - 2) THEN
               IZ = 0
               CALL FREFLD (IZ, IZ, 'TARGET ERROR MAGNITUDE,'//
     &            ' MIN FACTOR, 1 FACTOR, MAX FACTOR: ', MCOM,
     &            IOSTAT, JCOM, KIN, CIN, IIN, RIN)
               ICOM = 1
            ENDIF
            TERR = RIN (ICOM)
            ICOM = ICOM + 1
            EMINS = RIN (ICOM) * TERR
            ICOM = ICOM + 1
            E1S   = RIN (ICOM) * TERR
            ICOM = ICOM + 1
            EMAXS = RIN (ICOM) * TERR
            ICOM = ICOM + 1
         ENDIF
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 8/2/90
CC* MODIFICATION: ADDED INQUIRY FOR USER SPECIFIED CONTROL OF PAVING
CC**              REMESH PARAMETERS
C
C
C  GET THE MAXIMUM AND MINIMUM EXPANSION FACTORS
C
         IF (ICOM .GT. JCOM - 1) THEN
            IZ = 0
            CALL FREFLD (IZ, IZ, 'MAX EXPANSION FACTOR,'//
     &         ' MAX REDUCTION FACTOR: ', MCOM,
     &         IOSTAT, JCOM, KIN, CIN, IIN, RIN)
            ICOM = 1
         ENDIF
         EMAX = RIN (ICOM)
         ICOM = ICOM + 1
         EMIN = RIN (ICOM)
         ICOM = ICOM + 1
C
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/31/90
CC* MODIFICATION: ADDED MINIMUM ELEMENT SIZE FOR SINGULARITIES DURING
CC**              ERROR ESTIMATION.
C
C  CHECK TO SEE IF A MINIMUM ELEMENT SIZE IS NEEDED
C
         CALL INTRUP ('MINIMUM ELEMENT SIZE REQUIRED', MINSIZ, MCOM,
     &      ICOM, JCOM, CIN, IIN, RIN, KIN)
         IF (MINSIZ) THEN
            IF (ICOM .GT. JCOM) THEN
               IZ = 0
               CALL FREFLD (IZ, IZ, 'MINIMUM ELEMENT SIZE: ', MCOM,
     &            IOSTAT, JCOM, KIN, CIN, IIN, RIN)
               ICOM = 1
            ENDIF
            SIZMIN = RIN (ICOM)
            ICOM = ICOM + 1
         ELSE
            SIZMIN = 0.
         ENDIF
C
         FOUND1 = .FALSE.
         FOUND2 = .FALSE.
C
C  NOW GET THE TRUTH TABLE READ IN AND THEN READ EACH TIME STEP
C
         READ (IUNIT, ERR = 270) ((ISEVOK(I,J), I = 1, NVAREL),
     &      J = 1, NUMMAT)
         READ (IUNIT, ERR = 270) TIME, IFLAG
         READ (IUNIT, ERR = 270) (XDUM, I = 1, NVARHI)
C
         WRITE (IIUNIT, ERR = 270) ((ISEVOK(I,J), I = 1, 0),
     &      J = 1, NUMMAT)
         WRITE (IIUNIT, ERR = 270) TIME, IFLAG
         WRITE (IIUNIT, ERR = 270) (XDUM, I = 1, 0)
         IF (IFLAG .EQ. 0) THEN
            READ (IUNIT, ERR = 270) (XDUM, I = 1, NVARHI)
            WRITE (IIUNIT, ERR = 270) (XDUM, I = 1, 0)
C
C  READ IN THE NODAL VARIABLES
C
            DO 180 I = 1, NVARNP
               IF (IVAR1 .EQ. I) THEN
                  READ (IUNIT, ERR = 270) (BMESUR (J), J = 1, NNN)
                  FOUND1 = .TRUE.
               ELSEIF (IVAR2 .EQ. I) THEN
                  READ (IUNIT, ERR = 270) (DMESUR (J), J = 1, NNN)
                  FOUND2 = .TRUE.
               ELSE
                  READ (IUNIT, ERR = 270) (XDUM, J = 1, NNN)
               ENDIF
               IF ( ((FOUND1) .AND. (.NOT. SECOND)) .OR.
     &            ((FOUND1) .AND. (FOUND2)) ) GOTO 210
  180       CONTINUE
C
C  READ IN THE ELEMENT VARIABLES
C
            DO 200 IMAT = 1, NUMMAT
               DO 190 IEVAR = 1, NVAREL
                  IF (ISEVOK (IEVAR, IMAT) .NE. 0) THEN
                     IF (IVAR1 .EQ. IEVAR + NVARNP) THEN
                        READ (IUNIT, ERR = 270) (AMESUR (J),
     &                     J = MATMAP (2, IMAT), MATMAP (3, IMAT))
                        FOUND1 = .TRUE.
                     ELSEIF (IVAR2 .EQ. IEVAR + NVARNP) THEN
                        READ (IUNIT, ERR = 270) (CMESUR (J),
     &                     J = MATMAP (2, IMAT), MATMAP (3, IMAT))
                        FOUND2 = .TRUE.
                     ELSE
                        READ (IUNIT, ERR = 270) (XDUM,
     &                     J = MATMAP (2, IMAT), MATMAP (3, IMAT))
                     ENDIF
                     IF ( ((FOUND1) .AND. (.NOT. SECOND)) .OR.
     &                  ((FOUND1) .AND. (FOUND2)) ) GOTO 210
                  ENDIF
  190          CONTINUE
  200       CONTINUE
            CALL MESAGE ('** VARIABLE(S) NOT FOUND IN EXODUS **')
            GOTO 270
         ELSE
            CALL MESAGE ('** PROBLEMS WITH NO RESULTS IN EXODUS **')
            GOTO 270
         ENDIF
C
C  CONVERT THE ELEMENT VARIABLE(S) TO NODAL VARIABLE(S)
C
  210    CONTINUE
         CALL MDRSRV ('KNODE', K, NPNODE)
         IF (IVAR1 .GT. NVARNP)
     &      CALL AMXBM (NPNODE, NPELEM, NXK, AMESUR, BMESUR, IA(K))
         IF ((SECOND) .AND. (IVAR2 .GT. NVARNP))
     &      CALL AMXBM (NPNODE, NPELEM, NXK, CMESUR, DMESUR, IA(K))
         CALL MDDEL ('KNODE')
         CALL MDSTAT (NERR, MUSED)
         IF (NERR.GT.0) THEN
            CALL MDEROR (6)
            STOP ' '
         END IF
C
C  SET UP THE MAXIMUM AND MINIMUM REDUCTION FACTORS AS WELL AS THE TARGET
C  ERROR PARAMETERS
C
C  NORMALIZE AND COMBINE THE TWO VARIABLES (IF SUCH IS REQUIRED)
C
         CALL MIXVAR (NPNODE, BMESUR, DMESUR, SECOND, TERROR,
     &      EMIN, EMAX, EMINS, EMAXS, E1S)
         WRITE (IIUNIT, ERR = 270) (BMESUR (J), J = 1, NNN)
         DO 220 I = 1, NNN
            DMESUR(I) = 1.2 - BMESUR(I)*1.0571429
  220    CONTINUE
         WRITE (IIUNIT, ERR = 270) (DMESUR (J), J = 1, NNN)
         CLOSE (IIUNIT)
C
C  SET UP THE SPACE MAPPING OF THE ELEMENTS
C
         CALL ELARAY (XN, YN, NXK, MATMAP, LINKEG, LISTEG,
     &      MLINK, NPREGN, NPNOLD, NPEOLD, 9, REXMIN, REXMAX,
     &      REYMIN, REYMAX, IDIVIS)
      ENDIF

C
C  SET UP THE MAPGXD ARRAY
C
      DO 230 I = 1, KKK
         MAPGXD (MAPDXG (I)) = I
  230 CONTINUE
C
C  SET UP THE ELEMENT MATERIAL ARRAY
C
      DO 250 I = 1, NUMMAT
         DO 240 J = MATMAP (2, I), MATMAP (3, I)
            MAT (J) = MATMAP (1, I)
  240    CONTINUE
  250 CONTINUE
C
C  SUCCESSFULL READ COMPLETED
C
      CALL MESAGE (' ')
      CALL MESAGE ('**************************************************')
      CALL MESAGE ('**                                              **')
      CALL MESAGE ('**         GENESIS FILE SUCCESSFULLY READ       **')
      LWID = 0
      DO 260 K = 1, KKK
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
  260 CONTINUE
      WRITE(*, 10010) LWID
      WRITE(*, 10020) NNN, KKK, NUMMAT
      CALL MESAGE ('**                                              **')
      CALL MESAGE ('**************************************************')
C
      ERR = .FALSE.
      END = .FALSE.
      RETURN
C
C  ERR DURING READ
C
  270 CONTINUE
      CALL MESAGE ('ERR DURING READ OF GENESIS FILE')
      END = .FALSE.
      RETURN
C
10000 FORMAT(10X,I5,') ',A8)
10010 FORMAT(' **   LARGEST NODE DIFFERENCE PER ELEMENT:', I5, '  **')
10020 FORMAT(' ** NODES:', I5, '; ELEMENTS:', I5, '; MATERIALS:', I5,
     &   ' **')
C
      END
