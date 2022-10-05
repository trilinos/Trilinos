C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      SUBROUTINE HEXSRC(
     *  NDIM,     NPTS,     NPSRF,    NFSRF,    NISR,
     *  NRSR,     NRSS,     XYZE,     XYZP,     LS,
     *  ISRCHR,   RSRCHR,   IPT,      IELT,     IERR    )

C-----------------------------------------------------------------------

C DESCRIPTION:

C THIS SUBROUTINE CALCULATES THE CLOSEST POINT PROBLEM
C BETWEEN PAIRS OF POINTS AND SURFACES.

C-----------------------------------------------------------------------

C FORMAL PARAMETERS

C MEMORY      : P=PERMANENT, S=SCRATCH
C NAME        : IMPLICIT A-H,O-Z REAL, I-N INTEGER
C TYPE        : INPUT_STATUS/OUTPUT_STATUS (I=INPUT,O=OUTPUT,P=PASSED,
C               U=UNMODIFIED,-=UNDEFINED)
C DESCRIPTION : DESCRIPTION OF VARIABLE

C-----------------------------------------------------------------------

C CALLING ARGUMENTS

C MEMORY NAME     TYPE   DESCRIPTION
C ---    ----     ---    -----------
C  P     NDIM     I/U    DIMENSION OF PROBLEM=3
C  P     NPTS     I/U    NUMBER OF POINTS TO BE SEARCHED
C  P     NPSRF    I/U    NUMBER OF POINTS THAT DEFINE THE SURFACE
C  P     NFSRF    I/U    NUMBER OF SURFACES
C  P     NISR     I/U    NUMBER OF INTEGER SEARCH RESULTS (>=1)
C  P     NRSR     I/U    NUMBER OF REAL SEARCH RESULTS (>=4)
C  P     NRSS     I/U    NUMBER OF REAL SEARCH SCRATCH MEMORY (=10)
C  P     XYZE     I/U    XYZ COORDS OF POINTS DEFINING ELEMENT
C  P     XYZP     I/U    XYZ COORDS OF POINTS TO BE SEARCHED
C  P     LS       I/U    CONNECTIVITY OF ELEMENTS (8*NFSRF),
C                        NUMBERS REFER TO LOCATIONS IN XYZE ARRAY
C  P     ISRCHR   I/O    INTEGER SEARCH RESULTS
C  P     RSRCHR   I/O    REAL SEARCH RESULTS
C  P     IPT      I/U    POINT PAIRED WITH SURFACE LISTED IN IELT
C  P     IELT     I/U    SURFACE PAIRED WITH POINT LISTED IN IPT

C-----------------------------------------------------------------------

      include 'toldat.blk'
      include 'tapes.blk'

C INPUT/OUTPUT ARRAYS
      DIMENSION
     *  XYZP(NPTS,NDIM)     ,XYZE(NPSRF,NDIM)  ,LS(8,NFSRF)    ,
     *  ISRCHR(NISR,NPTS)   ,RSRCHR(NRSR,NPTS)
      DIMENSION XX(27), YY(27), ZZ(27)

      IF( NISR .LT. 1 .OR. NRSR .LT. 3 .OR. NRSS .LT. 10 )THEN
        IERR = 1
        RETURN
      ENDIF

C check for Mesh-B point coincident with node of element in Mesh-A

      SIDE1 = (XYZE(LS(1,IELT),1)-XYZE(LS(2,IELT),1))**2
     &      + (XYZE(LS(1,IELT),2)-XYZE(LS(2,IELT),2))**2
     &      + (XYZE(LS(1,IELT),3)-XYZE(LS(2,IELT),3))**2
      SIDE2 = (XYZE(LS(2,IELT),1)-XYZE(LS(3,IELT),1))**2
     &      + (XYZE(LS(2,IELT),2)-XYZE(LS(3,IELT),2))**2
     &      + (XYZE(LS(2,IELT),3)-XYZE(LS(3,IELT),3))**2
      SIDE3 = (XYZE(LS(3,IELT),1)-XYZE(LS(4,IELT),1))**2
     &      + (XYZE(LS(3,IELT),2)-XYZE(LS(4,IELT),2))**2
     &      + (XYZE(LS(3,IELT),3)-XYZE(LS(4,IELT),3))**2
      SIDE4 = (XYZE(LS(4,IELT),1)-XYZE(LS(1,IELT),1))**2
     &      + (XYZE(LS(4,IELT),2)-XYZE(LS(1,IELT),2))**2
     &      + (XYZE(LS(4,IELT),3)-XYZE(LS(1,IELT),3))**2
      SIDE5 = (XYZE(LS(5,IELT),1)-XYZE(LS(6,IELT),1))**2
     &      + (XYZE(LS(5,IELT),2)-XYZE(LS(6,IELT),2))**2
     &      + (XYZE(LS(5,IELT),3)-XYZE(LS(6,IELT),3))**2
      SIDE6 = (XYZE(LS(6,IELT),1)-XYZE(LS(7,IELT),1))**2
     &      + (XYZE(LS(6,IELT),2)-XYZE(LS(7,IELT),2))**2
     &      + (XYZE(LS(6,IELT),3)-XYZE(LS(7,IELT),3))**2
      SIDE7 = (XYZE(LS(7,IELT),1)-XYZE(LS(8,IELT),1))**2
     &      + (XYZE(LS(7,IELT),2)-XYZE(LS(8,IELT),2))**2
     &      + (XYZE(LS(7,IELT),3)-XYZE(LS(8,IELT),3))**2
      SIDE8 = (XYZE(LS(8,IELT),1)-XYZE(LS(5,IELT),1))**2
     &      + (XYZE(LS(8,IELT),2)-XYZE(LS(5,IELT),2))**2
     &      + (XYZE(LS(8,IELT),3)-XYZE(LS(5,IELT),3))**2
      SIDE9 = (XYZE(LS(1,IELT),1)-XYZE(LS(5,IELT),1))**2
     &      + (XYZE(LS(1,IELT),2)-XYZE(LS(5,IELT),2))**2
     &      + (XYZE(LS(1,IELT),3)-XYZE(LS(5,IELT),3))**2
      SIDE10 = (XYZE(LS(2,IELT),1)-XYZE(LS(6,IELT),1))**2
     &       + (XYZE(LS(2,IELT),2)-XYZE(LS(6,IELT),2))**2
     &       + (XYZE(LS(2,IELT),3)-XYZE(LS(6,IELT),3))**2
      SIDE11 = (XYZE(LS(3,IELT),1)-XYZE(LS(7,IELT),1))**2
     &       + (XYZE(LS(3,IELT),2)-XYZE(LS(7,IELT),2))**2
     &       + (XYZE(LS(3,IELT),3)-XYZE(LS(7,IELT),3))**2
      SIDE12 = (XYZE(LS(4,IELT),1)-XYZE(LS(8,IELT),1))**2
     &       + (XYZE(LS(4,IELT),2)-XYZE(LS(8,IELT),2))**2
     &       + (XYZE(LS(4,IELT),3)-XYZE(LS(8,IELT),3))**2
      SIDMAX = MAX(SIDE1,SIDE2,SIDE3,SIDE4,SIDE5,SIDE6,SIDE7,SIDE8
     &            ,SIDE9,SIDE10,SIDE11,SIDE12)
      SIDMIN = MIN(SIDE1,SIDE2,SIDE3,SIDE4,SIDE5,SIDE6,SIDE7,SIDE8
     &            ,SIDE9,SIDE10,SIDE11,SIDE12)

      COTEST = EPS*EPS*EPS*SIDMAX
      PRITST = COTEST

C test value for coincident nodes - used for mapping and checking
C                                   for prism element

C first check if recipient node is coincident with donor node
C if so no need to do anything else

      DO 110 I = 1, 8
        A = XYZE(LS(I,IELT),1) - XYZP(IPT,1)
        B = XYZE(LS(I,IELT),2) - XYZP(IPT,2)
        C = XYZE(LS(I,IELT),3) - XYZP(IPT,3)
        DIST = A**2 + B**2 + C**2
        IF (DIST .LT. COTEST)THEN

C coincident node, so fill search results arrays
C no need to check for better search result

          INODE = I
          ISRCHR(1,IPT) = IELT
          CALL NODE (10,INODE,RSRCHR(1,IPT),RSRCHR(2,IPT),
     &      RSRCHR(3,IPT))
          GO TO 100
        END IF
 110  CONTINUE

C check  is hex degenerate

      IF (SIDMIN .LT. COTEST)THEN

C hex is degenerate,  check is it a prism
C if prism, build prism jacobian

        IF (SIDE1 .LT. PRITST .AND. SIDE3 .LT. PRITST)THEN
          XX(1) = XYZE(LS(3,IELT),1)
          XX(2) = XYZE(LS(7,IELT),1)
          XX(3) = XYZE(LS(8,IELT),1)
          XX(4) = XYZE(LS(1,IELT),1)
          XX(5) = XYZE(LS(6,IELT),1)
          XX(6) = XYZE(LS(5,IELT),1)
          YY(1) = XYZE(LS(3,IELT),2)
          YY(2) = XYZE(LS(7,IELT),2)
          YY(3) = XYZE(LS(8,IELT),2)
          YY(4) = XYZE(LS(1,IELT),2)
          YY(5) = XYZE(LS(6,IELT),2)
          YY(6) = XYZE(LS(5,IELT),2)
          ZZ(1) = XYZE(LS(3,IELT),3)
          ZZ(2) = XYZE(LS(7,IELT),3)
          ZZ(3) = XYZE(LS(8,IELT),3)
          ZZ(4) = XYZE(LS(1,IELT),3)
          ZZ(5) = XYZE(LS(6,IELT),3)
          ZZ(6) = XYZE(LS(5,IELT),3)
          GO TO 200
        ELSE IF (SIDE1 .LT. PRITST .AND. SIDE5 .LT. PRITST)THEN
          XX(1) = XYZE(LS(1,IELT),1)
          XX(2) = XYZE(LS(3,IELT),1)
          XX(3) = XYZE(LS(4,IELT),1)
          XX(4) = XYZE(LS(5,IELT),1)
          XX(5) = XYZE(LS(7,IELT),1)
          XX(6) = XYZE(LS(8,IELT),1)
          YY(1) = XYZE(LS(1,IELT),2)
          YY(2) = XYZE(LS(3,IELT),2)
          YY(3) = XYZE(LS(4,IELT),2)
          YY(4) = XYZE(LS(5,IELT),2)
          YY(5) = XYZE(LS(7,IELT),2)
          YY(6) = XYZE(LS(8,IELT),2)
          YY(1) = XYZE(LS(1,IELT),2)
          YY(2) = XYZE(LS(3,IELT),2)
          YY(3) = XYZE(LS(4,IELT),2)
          YY(4) = XYZE(LS(5,IELT),2)
          YY(5) = XYZE(LS(7,IELT),2)
          YY(6) = XYZE(LS(8,IELT),2)
          GO TO 200
        ELSE IF (SIDE2 .LT. PRITST .AND. SIDE4 .LT. PRITST)THEN
          XX(1) = XYZE(LS(2,IELT),1)
          XX(2) = XYZE(LS(6,IELT),1)
          XX(3) = XYZE(LS(7,IELT),1)
          XX(4) = XYZE(LS(1,IELT),1)
          XX(5) = XYZE(LS(5,IELT),1)
          XX(6) = XYZE(LS(8,IELT),1)
          YY(1) = XYZE(LS(2,IELT),2)
          YY(2) = XYZE(LS(6,IELT),2)
          YY(3) = XYZE(LS(7,IELT),2)
          YY(4) = XYZE(LS(1,IELT),2)
          YY(5) = XYZE(LS(5,IELT),2)
          YY(6) = XYZE(LS(8,IELT),2)
          ZZ(1) = XYZE(LS(2,IELT),3)
          ZZ(2) = XYZE(LS(6,IELT),3)
          ZZ(3) = XYZE(LS(7,IELT),3)
          ZZ(4) = XYZE(LS(1,IELT),3)
          ZZ(5) = XYZE(LS(5,IELT),3)
          ZZ(6) = XYZE(LS(8,IELT),3)
          GO TO 200
        ELSE IF (SIDE2 .LT. PRITST .AND. SIDE6 .LT. PRITST)THEN
          XX(1) = XYZE(LS(2,IELT),1)
          XX(2) = XYZE(LS(4,IELT),1)
          XX(3) = XYZE(LS(1,IELT),1)
          XX(4) = XYZE(LS(6,IELT),1)
          XX(5) = XYZE(LS(8,IELT),1)
          XX(6) = XYZE(LS(5,IELT),1)
          YY(1) = XYZE(LS(2,IELT),2)
          YY(2) = XYZE(LS(4,IELT),2)
          YY(3) = XYZE(LS(1,IELT),2)
          YY(4) = XYZE(LS(6,IELT),2)
          YY(5) = XYZE(LS(8,IELT),2)
          YY(6) = XYZE(LS(5,IELT),2)
          ZZ(1) = XYZE(LS(2,IELT),3)
          ZZ(2) = XYZE(LS(4,IELT),3)
          ZZ(3) = XYZE(LS(1,IELT),3)
          ZZ(4) = XYZE(LS(6,IELT),3)
          ZZ(5) = XYZE(LS(8,IELT),3)
          ZZ(6) = XYZE(LS(5,IELT),3)
          GO TO 200
        ELSE IF (SIDE3 .LT. PRITST .AND. SIDE7 .LT. PRITST)THEN
          XX(1) = XYZE(LS(3,IELT),1)
          XX(2) = XYZE(LS(1,IELT),1)
          XX(3) = XYZE(LS(2,IELT),1)
          XX(4) = XYZE(LS(7,IELT),1)
          XX(5) = XYZE(LS(5,IELT),1)
          XX(6) = XYZE(LS(6,IELT),1)
          YY(1) = XYZE(LS(3,IELT),2)
          YY(2) = XYZE(LS(1,IELT),2)
          YY(3) = XYZE(LS(2,IELT),2)
          YY(4) = XYZE(LS(7,IELT),2)
          YY(5) = XYZE(LS(5,IELT),2)
          YY(6) = XYZE(LS(6,IELT),2)
          ZZ(1) = XYZE(LS(3,IELT),3)
          ZZ(2) = XYZE(LS(1,IELT),3)
          ZZ(3) = XYZE(LS(2,IELT),3)
          ZZ(4) = XYZE(LS(7,IELT),3)
          ZZ(5) = XYZE(LS(5,IELT),3)
          ZZ(6) = XYZE(LS(6,IELT),3)
          GO TO 200
        ELSE IF (SIDE4 .LT. PRITST .AND. SIDE8 .LT. PRITST)THEN
          XX(1) = XYZE(LS(1,IELT),1)
          XX(2) = XYZE(LS(2,IELT),1)
          XX(3) = XYZE(LS(3,IELT),1)
          XX(4) = XYZE(LS(5,IELT),1)
          XX(5) = XYZE(LS(6,IELT),1)
          XX(6) = XYZE(LS(7,IELT),1)
          YY(1) = XYZE(LS(1,IELT),2)
          YY(2) = XYZE(LS(2,IELT),2)
          YY(3) = XYZE(LS(3,IELT),2)
          YY(4) = XYZE(LS(5,IELT),2)
          YY(5) = XYZE(LS(6,IELT),2)
          YY(6) = XYZE(LS(7,IELT),2)
          ZZ(1) = XYZE(LS(1,IELT),3)
          ZZ(2) = XYZE(LS(2,IELT),3)
          ZZ(3) = XYZE(LS(3,IELT),3)
          ZZ(4) = XYZE(LS(5,IELT),3)
          ZZ(5) = XYZE(LS(6,IELT),3)
          ZZ(6) = XYZE(LS(7,IELT),3)
          GO TO 200
        ELSE IF (SIDE5 .LT. PRITST .AND. SIDE7 .LT. PRITST)THEN
          XX(1) = XYZE(LS(5,IELT),1)
          XX(2) = XYZE(LS(2,IELT),1)
          XX(3) = XYZE(LS(1,IELT),1)
          XX(4) = XYZE(LS(7,IELT),1)
          XX(5) = XYZE(LS(3,IELT),1)
          XX(6) = XYZE(LS(4,IELT),1)
          YY(1) = XYZE(LS(5,IELT),2)
          YY(2) = XYZE(LS(2,IELT),2)
          YY(3) = XYZE(LS(1,IELT),2)
          YY(4) = XYZE(LS(7,IELT),2)
          YY(5) = XYZE(LS(3,IELT),2)
          YY(6) = XYZE(LS(4,IELT),2)
          ZZ(1) = XYZE(LS(5,IELT),3)
          ZZ(2) = XYZE(LS(2,IELT),3)
          ZZ(3) = XYZE(LS(1,IELT),3)
          ZZ(4) = XYZE(LS(7,IELT),3)
          ZZ(5) = XYZE(LS(3,IELT),3)
          ZZ(6) = XYZE(LS(4,IELT),3)
          GO TO 200
        ELSE IF (SIDE6 .LT. PRITST .AND. SIDE8 .LT. PRITST)THEN
          XX(1) = XYZE(LS(5,IELT),1)
          XX(2) = XYZE(LS(1,IELT),1)
          XX(3) = XYZE(LS(4,IELT),1)
          XX(4) = XYZE(LS(7,IELT),1)
          XX(5) = XYZE(LS(2,IELT),1)
          XX(6) = XYZE(LS(3,IELT),1)
          YY(1) = XYZE(LS(5,IELT),2)
          YY(2) = XYZE(LS(1,IELT),2)
          YY(3) = XYZE(LS(4,IELT),2)
          YY(4) = XYZE(LS(6,IELT),2)
          YY(5) = XYZE(LS(2,IELT),2)
          YY(6) = XYZE(LS(3,IELT),2)
          ZZ(1) = XYZE(LS(5,IELT),3)
          ZZ(2) = XYZE(LS(1,IELT),3)
          ZZ(3) = XYZE(LS(4,IELT),3)
          ZZ(4) = XYZE(LS(6,IELT),3)
          ZZ(5) = XYZE(LS(2,IELT),3)
          ZZ(6) = XYZE(LS(3,IELT),3)
          GO TO 200
        ELSE IF (SIDE9 .LT. PRITST .AND. SIDE10 .LT. PRITST)THEN
          XX(1) = XYZE(LS(1,IELT),1)
          XX(2) = XYZE(LS(4,IELT),1)
          XX(3) = XYZE(LS(8,IELT),1)
          XX(4) = XYZE(LS(2,IELT),1)
          XX(5) = XYZE(LS(3,IELT),1)
          XX(6) = XYZE(LS(7,IELT),1)
          YY(1) = XYZE(LS(1,IELT),2)
          YY(2) = XYZE(LS(4,IELT),2)
          YY(3) = XYZE(LS(8,IELT),2)
          YY(4) = XYZE(LS(2,IELT),2)
          YY(5) = XYZE(LS(3,IELT),2)
          YY(6) = XYZE(LS(7,IELT),2)
          ZZ(1) = XYZE(LS(1,IELT),3)
          ZZ(2) = XYZE(LS(4,IELT),3)
          ZZ(3) = XYZE(LS(8,IELT),3)
          ZZ(4) = XYZE(LS(2,IELT),3)
          ZZ(5) = XYZE(LS(3,IELT),3)
          ZZ(6) = XYZE(LS(7,IELT),3)
          GO TO 200
        ELSE IF (SIDE10 .LT. PRITST .AND. SIDE11 .LT. PRITST)THEN
          XX(1) = XYZE(LS(2,IELT),1)
          XX(2) = XYZE(LS(1,IELT),1)
          XX(3) = XYZE(LS(5,IELT),1)
          XX(4) = XYZE(LS(3,IELT),1)
          XX(5) = XYZE(LS(4,IELT),1)
          XX(6) = XYZE(LS(8,IELT),1)
          YY(1) = XYZE(LS(2,IELT),2)
          YY(2) = XYZE(LS(1,IELT),2)
          YY(3) = XYZE(LS(5,IELT),2)
          YY(4) = XYZE(LS(3,IELT),2)
          YY(5) = XYZE(LS(4,IELT),2)
          YY(6) = XYZE(LS(8,IELT),2)
          ZZ(1) = XYZE(LS(2,IELT),3)
          ZZ(2) = XYZE(LS(1,IELT),3)
          ZZ(3) = XYZE(LS(5,IELT),3)
          ZZ(4) = XYZE(LS(3,IELT),3)
          ZZ(5) = XYZE(LS(4,IELT),3)
          ZZ(6) = XYZE(LS(8,IELT),3)
          GO TO 200
        ELSE IF (SIDE11 .LT. PRITST .AND. SIDE12 .LT. PRITST)THEN
          XX(1) = XYZE(LS(3,IELT),1)
          XX(2) = XYZE(LS(2,IELT),1)
          XX(3) = XYZE(LS(6,IELT),1)
          XX(4) = XYZE(LS(4,IELT),1)
          XX(5) = XYZE(LS(1,IELT),1)
          XX(6) = XYZE(LS(5,IELT),1)
          YY(1) = XYZE(LS(3,IELT),2)
          YY(2) = XYZE(LS(2,IELT),2)
          YY(3) = XYZE(LS(6,IELT),2)
          YY(4) = XYZE(LS(4,IELT),2)
          YY(5) = XYZE(LS(1,IELT),2)
          YY(6) = XYZE(LS(5,IELT),2)
          ZZ(1) = XYZE(LS(3,IELT),3)
          ZZ(2) = XYZE(LS(2,IELT),3)
          ZZ(3) = XYZE(LS(6,IELT),3)
          ZZ(4) = XYZE(LS(4,IELT),3)
          ZZ(5) = XYZE(LS(1,IELT),3)
          ZZ(6) = XYZE(LS(5,IELT),3)
          GO TO 200
        ELSE IF (SIDE12 .LT. PRITST .AND. SIDE9 .LT. PRITST)THEN
          XX(1) = XYZE(LS(4,IELT),1)
          XX(2) = XYZE(LS(3,IELT),1)
          XX(3) = XYZE(LS(7,IELT),1)
          XX(4) = XYZE(LS(1,IELT),1)
          XX(5) = XYZE(LS(2,IELT),1)
          XX(6) = XYZE(LS(6,IELT),1)
          YY(1) = XYZE(LS(4,IELT),2)
          YY(2) = XYZE(LS(3,IELT),2)
          YY(3) = XYZE(LS(7,IELT),2)
          YY(4) = XYZE(LS(1,IELT),2)
          YY(5) = XYZE(LS(2,IELT),2)
          YY(6) = XYZE(LS(6,IELT),2)
          ZZ(1) = XYZE(LS(4,IELT),3)
          ZZ(2) = XYZE(LS(3,IELT),3)
          ZZ(3) = XYZE(LS(7,IELT),3)
          ZZ(4) = XYZE(LS(1,IELT),3)
          ZZ(5) = XYZE(LS(2,IELT),3)
          ZZ(6) = XYZE(LS(6,IELT),3)
          GO TO 200
        ELSE
          CALL ERROR ('SRCHH',
     &      'DEGENERATE HEX',
     &      'MESH-A ELEMENT',IELT,
     &      'IS NOT A PRISM',0,' ',' ',0)
          GO TO 100
        END IF
 200    CONTINUE

C process  prism element

        SG = 0.
        TG = 0.
        RG = 0.
        ITER = 0

 210    CONTINUE
        CALL JACOBN (8,XX,YY,ZZ,SG,TG,RG,A11,A12,A13,A21,A22,A23,
     &    A31,A32,A33,F1,F2,F3)
        DETA = A11*(A22*A33 - A23*A32) - A12*(A21*A33 - A31*A23)
     &    + A13*(A21*A32 - A31*A22)
        IF (ABS(DETA) .LT. 1.E-25)THEN
          CALL ERROR ('SRCHH',
     &      'ZERO JACOBIAN FOUND DURING NEWTON ITERATION',
     &      'MESH-A ELEMENT',IELT,
     &      'TRYING TO PROCESS AS A PRISM',
     &      0,' ',' ',0)
C diverging Newton iteration
          GOTO 100
        END IF

        AI11 =  (A22*A33 - A23*A32)/DETA
        AI12 = -(A12*A33 - A32*A13)/DETA
        AI13 =  (A23*A12 - A13*A22)/DETA
        AI21 = -(A21*A33 - A31*A23)/DETA
        AI22 =  (A11*A33 - A13*A31)/DETA
        AI23 = -(A11*A23 - A21*A13)/DETA
        AI31 =  (A21*A32 - A31*A22)/DETA
        AI32 = -(A11*A32 - A31*A12)/DETA
        AI33 =  (A11*A22 - A12*A21)/DETA

        FS = F1 - XYZP(IPT,1)
        FT = F2 - XYZP(IPT,2)
        FR = F3 - XYZP(IPT,3)
        SNEW = SG - (AI11*FS + AI12*FT + AI13*FR)
        TNEW = TG - (AI21*FS + AI22*FT + AI23*FR)
        RNEW = RG - (AI31*FS + AI32*FT + AI33*FR)

        ITER = ITER + 1
        DS = ABS(SNEW-SG)
        DT = ABS(TNEW-TG)
        DR = ABS(RNEW-RG)
        IF (DS .LT. TOL .AND. DT .LT. TOL .AND. DR .LT. TOL) GO TO 300
        SG = SNEW
        TG = TNEW
        RG = RNEW
C check for diverging Newton iteration
        IF ( (ITER .GE. ITERMX) .OR.
     *     ( ABS( SG) .GT. 1E+06) .OR.
     *     ( ABS( TG) .GT. 1E+06) .OR.
     *     ( ABS( RG) .GT. 1E+06)) GOTO 100
        GO TO 210
      ELSE

C hex element
C compute isoparametric coordinates for hex - Use Newton's method

        SG = 0.
        TG = 0.
        RG = 0.
        ITER = 0

C Build Jacobian and invert

        DO 120 I = 1, 8
          XX(I) = XYZE(LS(I,IELT),1)
          YY(I) = XYZE(LS(I,IELT),2)
          ZZ(I) = XYZE(LS(I,IELT),3)
 120    CONTINUE
 130    CONTINUE
        CALL JACOBN (10,XX,YY,ZZ,SG,TG,RG,A11,A12,A13,A21,A22,A23,
     &    A31,A32,A33,F1,F2,F3)
        DETA = A11*(A22*A33 - A23*A32) - A12*(A21*A33 - A31*A23)
     &    + A13*(A21*A32 - A31*A22)
        IF (ABS(DETA) .LT. 1.E-25)THEN
          CALL ERROR ('SRCHH',
     &      'ZERO JACOBIAN FOUND DURING NEWTON ITERATION',
     &      'MESH-A ELEMENT',IELT,
     &      'TRYING TO PROCESS AS A HEX',
     &      0,' ',' ',0)
          GOTO 100
        END IF

        AI11 =  (A22*A33 - A23*A32)/DETA
        AI12 = -(A12*A33 - A32*A13)/DETA
        AI13 =  (A23*A12 - A13*A22)/DETA
        AI21 = -(A21*A33 - A31*A23)/DETA
        AI22 =  (A11*A33 - A13*A31)/DETA
        AI23 = -(A11*A23 - A21*A13)/DETA
        AI31 =  (A21*A32 - A31*A22)/DETA
        AI32 = -(A11*A32 - A31*A12)/DETA
        AI33 =  (A11*A22 - A12*A21)/DETA

        FS = F1 - XYZP(IPT,1)
        FT = F2 - XYZP(IPT,2)
        FR = F3 - XYZP(IPT,3)
        SNEW = SG - (AI11*FS + AI12*FT + AI13*FR)
        TNEW = TG - (AI21*FS + AI22*FT + AI23*FR)
        RNEW = RG - (AI31*FS + AI32*FT + AI33*FR)

        ITER = ITER + 1
        DS = ABS(SNEW-SG)
        DT = ABS(TNEW-TG)
        DR = ABS(RNEW-RG)
        IF (DS .LT. TOL .AND. DT .LT. TOL .AND. DR .LT. TOL) GO TO 300
        SG = SNEW
        TG = TNEW
        RG = RNEW
C check for diverging Newton iteration
        IF ( (ITER .GE. ITERMX) .OR.
     *     ( ABS( SG) .GT. 1E+06) .OR.
     *     ( ABS( TG) .GT. 1E+06) .OR.
     *     ( ABS( RG) .GT. 1E+06)) GOTO 100
        GO TO 130
      END IF

 300  CONTINUE

C Newton converged, load up search results arrays if appropriate

      IF (ABS(SNEW) .LT. STRLMT .AND. ABS(TNEW) .LT. STRLMT .AND.
     &  ABS(RNEW) .LT. STRLMT)THEN

C Search was adequate

        FTEST = MAX(ABS(RSRCHR(1,IPT)),ABS(RSRCHR(2,IPT)),
     &    ABS(RSRCHR(3,IPT)))
        FCOMP = MAX(ABS(SNEW),ABS(TNEW),ABS(RNEW))

        IF (FTEST .GT. FCOMP .OR. ISRCHR(1,IPT) .EQ. 0)THEN

C New search is better or first find, replace search results

          ISRCHR(1,IPT) = IELT
          RSRCHR(1,IPT) = SNEW
          RSRCHR(2,IPT) = TNEW
          RSRCHR(3,IPT) = RNEW
        END IF
      END IF
 100  CONTINUE
      RETURN
      END

