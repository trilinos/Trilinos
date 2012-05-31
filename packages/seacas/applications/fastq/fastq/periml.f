C $Id: periml.f,v 1.2 1998/07/14 18:19:33 gdsjaar Exp $
C $Log: periml.f,v $
C Revision 1.2  1998/07/14 18:19:33  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:13:21  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:13:19  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]PERIML.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE PERIML (NBNODE, MXND, NPER, ISTART, MLN, XN, YN, ZN,
     &   LXK, KXL, NXL, LXN, ANGLE, BNSIZE, LNODES, LPERIM, LLL, LLLOLD,
     &   XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, DEV1, KREG, ERR)
C***********************************************************************
C
C  SUBROUTINE PERIML = LINKS THE PERIMETER OF A REGION TOGETHER FOR
C                      THE FILL ROUTINES
C
C***********************************************************************
C
C  VARIABLES USED:
C     NPER  = NUMBER OF PERIMETER NODES
C     ERR   = .TRUE. IF ERRORS WERE ENCOUNTERED
C     XN    = GLOBAL X VALUES OF NODES
C     YN    = GLOBAL Y VALUES OF NODES
C     ZN    = GLOBAL Z VALUES OF NODES
C     LXK   = LINES PER ELEMENT
C     KXL   = ELEMENTS PER LINE
C     NXL   = NODES PER LINE
C     LXN   = LINES PER NODE
C  NOTE:
C     FOR *XN TABLES A NEGATIVE FLAG IN THE FOURTH COLUMN MEANS
C     GO TO THAT ROW FOR A CONTINUATION OF THE LIST.  IN THAT ROW
C     THE FIRST ELEMENT WILL BE NEGATED TO INDICATE THAT THIS IS
C     A CONTINUATION ROW.
C     A NEGATIVE FLAG IN THE SECOND COLUMN OF THE LXN ARRAY MEANS
C     THAT THIS NODE IS AN EXTERIOR BOUNDARY NODE.
C
C***********************************************************************
C
      DIMENSION XN (MXND), YN (MXND), ZN (MXND)
      DIMENSION LPERIM (NBNODE)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION ANGLE (MXND), BNSIZE (2, MXND), LNODES (MLN, MXND)
C
      CHARACTER*3 DEV1
C
      LOGICAL ERR
C
      ERR = .FALSE.
      LLLOLD = LLL
      IEND = NPER + ISTART - 1
C
      DO 100 I = ISTART, IEND
         NODE1 = LPERIM (I)
         IF (I .EQ. IEND) THEN
            NODE0 = LPERIM (I - 1)
            NODE2 = LPERIM (ISTART)
         ELSEIF (I .EQ. ISTART) THEN
            NODE0 = LPERIM (IEND)
            NODE2 = LPERIM (I + 1)
         ELSE
            NODE0 = LPERIM (I - 1)
            NODE2 = LPERIM (I + 1)
         ENDIF
         LLL = LLL+1
C
C  FILL UP THE NODES PER LINE ARRAY
C
         NXL (1, LLL) = NODE1
         NXL (2, LLL) = NODE2
C
C  FILL UP THE LINES PER NODE ARRAY
C
         LXN (1, NODE1) = LLL
         IF (I .EQ. ISTART) THEN
            LXN (2, NODE1) = - (IEND - ISTART + 1 + LLLOLD)
         ELSE
            LXN (2, NODE1) = 1 - LLL
         ENDIF
         LXN (3, NODE1) = 0
         LXN (4, NODE1) = 0
C
C  THE LNODES ARRAY IS DOCUMENTED IN THE ADDROW ROUTINE
C
         LNODES (1, NODE1) = 0
         LNODES (2, NODE1) = NODE0
         LNODES (3, NODE1) = NODE2
         LNODES (4, NODE1) = 1
         LNODES (5, NODE1) = LLL
         LNODES (6, NODE1) = 0
         LNODES (7, NODE1) = 0
         LNODES (8, NODE1) = 0
  100 CONTINUE
C
C  SET ALL THE INTERIOR ANGLES
C
      CALL LUPANG (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN, NPER,
     &   ANGLE, LNODES, LPERIM(ISTART), LLL, XMIN, XMAX, YMIN, YMAX,
     &   ZMIN, ZMAX, DEV1, KREG, ERR)
      IF (ERR) GOTO 120
C
C  SET UP THE CORRECT SIZE ARRAY
C
      N1S = LPERIM(ISTART)
      KOUNT = 0
  110 CONTINUE
      N0 = LNODES (2, N1S)
      N2 = LNODES (3, N1S)
      N1S = N2
      IF (N1S .EQ. LPERIM(ISTART)) GOTO 120
      KOUNT = KOUNT + 1
      IF (KOUNT .GT. NPER) THEN
         CALL MESAGE (' ** ERROR IN PERIML GETTING NODE SIZES ** ')
         ERR = .TRUE.
         GOTO 120
      ENDIF
C
      BNSIZE (1, N1S) =
     &   ( SQRT ( (XN (N1S) - XN (N0)) **2 + (YN (N1S) - YN (N0)) **2) +
     &   SQRT ( (XN (N1S) - XN (N2)) **2 + (YN (N1S) - YN (N2)) **2) )
     &   * .5
      BNSIZE (2, N1S) = 1.
      GOTO 110
C
  120 CONTINUE
C
      RETURN
C
      END
