C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      SUBROUTINE SRCHQ (
     *  NPSRF,   NFSRF,   LINKSRF, XYZSRF,
     *  NPTS,    XYZPTS,  TOLSRCH,
     *  NISR,    NRSR,    NISS,    NRSS,    ISRCHR,  RSRCHR, LBLK,
     *  LIST,    IND,     IRNK,    IRNK2,   INDX,    ILO,     IUP,
     *  IDP,     IDS,     XMIN,    XMAX,    ISCR,    RSCR,    IERR )

C-----------------------------------------------------------------------

C DESCRIPTION:

C THIS SUBROUTINE CALCULATES THE CLOSEST POINT PROBLEM
C BETWEEN NPTS POINTS AND NFSRF SURFACES AND RETURNS RESULTS OF
C SEARCH IN ISRCHR,RSRCHR

C USED HERE FOR FINDING LOCATION OF EITHER NODE OR ELEMENT CENTROID
C FROM MESH-B IN QUAD-4 ELEMENT OF MESH-A

C-----------------------------------------------------------------------

C FORMAL PARAMETERS

C MEMORY      : P=PERMANENT, S=SCRATCH
C NAME        : IMPLICIT A-H,O-Z REAL, I-N INTEGER
C TYPE        : INPUT_STATUS/OUTPUT_STATUS (I=INPUT,O=OUTPUT,P=PASSED,
C               U=UNMODIFIED,-=UNDEFINED)
C DESCRIPTION : DESCRIPTION OF VARIABLE

C-----------------------------------------------------------------------

C CALLING ARGUMENTS:

C MEMORY NAME     TYPE   DESCRIPTION
C ---    ----     ---    -----------
C  P     NPSRF    I/U    NUMBER OF POINTS THAT DEFINE THE SURFACE
C  P     NFSRF    I/U    NUMBER OF SURFACES
C  P     LINKSRF  I/U    CONNECTIVITY OF SURFACES OF SIZE (4*NFSRF),
C                        NUMBERS REFER TO LOCATIONS IN XYZSRF ARRAY
C  P     XYZSRF   I/U    XYZ COORDS OF POINTS DEFINING SURFACE
C  P     NPTS     I/U    NUMBER OF POINTS TO BE SEARCHED
C  P     XYZPTS   I/U    XYZ COORDS OF POINTS TO BE SEARCHED
C  P     TOLSRCH  I/U    PROXIMITY TOLERANCE FOR POINT-TO-SURFACE SEARCH
C  P     NISR     I/U    NUMBER OF INTEGER SEARCH RESULTS (>=1)
C  P     NRSR     I/U    NUMBER OF REAL SEARCH RESULTS (>=4)
C  P     NISS     I/U    NUMBER OF INTEGER SEARCH SCRATCH (=5)
C  P     NRSS     I/U    NUMBER OF REAL SEARCH SCRATCH (=10)
C  P     ISRCHR   -/O    INTEGER SEARCH RESULTS
C  P     RSRCHR   -/O    REAL SEARCH RESULTS
C  P     LBLK     I/U    VECTOR BLOCKING PARAMETER (RECOMMEND 512 ON CRAY)
C  S     LIST     -/-    LIST OF POTENTIAL CONTACTS FOR A SURFACE
C  S     IND      -/-    INDEX ARRAY INDICATING COORD. ORDER OF EACH POINT
C  S     IRNK     -/-    RANK ARRAY
C  S     IRNK2    -/-    RANK ARRAY (INDIRECT)
C  S     INDX     -/P    INTERMEDIATE LIST OF POTENTIAL PAIRS (CONSTRUCTED
C                          AFTER THE INTERSECTION OF TWO OF THREE LISTS)
C  S     ILO      -/-    SEARCH BOX MINIMUM INDEX
C  S     IUP      -/-    SEARCH BOX MAXIMUM INDEX
C  S     XMIN     -/-    SEARCH BOX MINIMUM DIMENSION
C  S     XMAX     -/-    SEARCH BOX MAXIMUM DIMENSION
C  S     ISCR     -/-    INTEGER SCRATCH MEMORY
C  S     RSCR     -/-    REAL SCRATCH MEMORY
C  P     IERR     -/O    ERROR FLAG

C-----------------------------------------------------------------------

      include 'amesh.blk'
      include 'tapes.blk'

C INPUT/OUTPUT ARRAYS

      DIMENSION
     *  LINKSRF(NELNDA,NFSRF),   XYZSRF(NPSRF,3),
     *  XYZPTS(NPTS,3),
     *  ISRCHR(NISR,NPTS),  RSRCHR(NRSR,NPTS)
      DIMENSION
     *  LIST(NPTS),         IND(NPTS,3),         IRNK(NPTS,3),
     *  IRNK2(NPTS,3,2),    INDX(NPTS)
      DIMENSION
     *  ILO(LBLK,3),        IUP(LBLK,3),         IDP(LBLK),
     *  IDS(LBLK)
      DIMENSION
     *  XMIN(LBLK,3),       XMAX(LBLK,3),        ISCR(NISS*LBLK),
     *  RSCR(NRSS*LBLK)

C ISRCHR and RSRCHR must be initialized to zero

      DO 1 I = 1, NPTS
        DO 2 J = 1, NISR
          ISRCHR(J,I) = 0
 2      CONTINUE
        DO 3 K = 1, NRSR
          RSRCHR(K,I) = 0.
 3      CONTINUE
 1    CONTINUE

      IF( NISR .LT. 1 .OR. NRSR .LT. 2 .OR. NISS .LT. 5 .OR.
     *    NRSS .LT. 10 )THEN
        IERR = 1
        RETURN
      ENDIF
C DIMENSION OF COORDINATES
      NDIM = 2

C CALL SORTING ROUTINE TO MAKE RANK ARRAYS
      CALL MKRNK( NPTS,NPTS,NDIM,XYZPTS,IND,IRNK,IRNK2 )

C LOOP OVER SURFACES AND SEARCH FOR POINTS WITHIN CAPTURE BOX
      DO 100 IFSRF = 1, NFSRF, LBLK
        NE = MIN(LBLK,NFSRF-IFSRF+1)
C CONSTRUCT THE BOUNDING BOX FOR NE SURFACES
        DO 110 J = 1, NE
          JFSRF = IFSRF + J - 1
          NI = LINKSRF(1,JFSRF)
          NJ = LINKSRF(2,JFSRF)
          NK = LINKSRF(3,JFSRF)
          NL = LINKSRF(4,JFSRF)

          XMINMS = MIN(XYZSRF(NI,1),XYZSRF(NJ,1),
     *                 XYZSRF(NK,1),XYZSRF(NL,1))
          XMAXMS = MAX(XYZSRF(NI,1),XYZSRF(NJ,1),
     *                 XYZSRF(NK,1),XYZSRF(NL,1))
          YMINMS = MIN(XYZSRF(NI,2),XYZSRF(NJ,2),
     *                 XYZSRF(NK,2),XYZSRF(NL,2))
          YMAXMS = MAX(XYZSRF(NI,2),XYZSRF(NJ,2),
     *                 XYZSRF(NK,2),XYZSRF(NL,2))
          ZMINMS = MIN(XYZSRF(NI,3),XYZSRF(NJ,3),
     *                 XYZSRF(NK,3),XYZSRF(NL,3))
          ZMAXMS = MAX(XYZSRF(NI,3),XYZSRF(NJ,3),
     *                 XYZSRF(NK,3),XYZSRF(NL,3))
          TOLER = TOLSRCH * (XMAXMS - XMINMS)
          XMIN(J,1) = XMINMS - TOLER
          XMAX(J,1) = XMAXMS + TOLER
          TOLER = TOLSRCH * (YMAXMS - YMINMS)
          XMIN(J,2) = YMINMS - TOLER
          XMAX(J,2) = YMAXMS + TOLER
          TOLER = TOLSRCH * (ZMAXMS - ZMINMS)
          XMIN(J,3) = ZMINMS - TOLER
          XMAX(J,3) = ZMAXMS + TOLER
  110   CONTINUE

        DO 120 IDIM = 1, NDIM
          CALL GETBND( LBLK, NE, XYZPTS(1,IDIM),  IND(1,IDIM),
     *                 NPTS,   XMIN(1,IDIM), XMAX(1,IDIM),
     *                 NPTS,    ILO(1,IDIM),  IUP(1,IDIM),
     *                 ISCR, RSCR )
  120   CONTINUE

        DO 130 J = 1, NE
          JFSRF = IFSRF + J - 1
          CALL MKLSTV(NPTS,IND,IRNK2,IUP,ILO,INDX,J,LIST,NLIST,
     *                LBLK,NDIM)

          DO 140 K = 1, NLIST
            lval = list(k)
            CALL QADSRC(
     *      NDIM,     NPTS,     NPSRF,    NFSRF,    NISR,
     *      NRSR,     NRSS,     XYZSRF,   XYZPTS,   LINKSRF,
     *      ISRCHR,   RSRCHR,   LVAL,     JFSRF,    IERR   )
  140     CONTINUE
  130   CONTINUE
  100 CONTINUE
      RETURN
      END

