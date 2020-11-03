C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      SUBROUTINE QADSRC(
     *  NDIM,     NPTS,     NPSRF,    NFSRF,    NISR,
     *  NRSR,     NRSS,     XYZE,     XYZP,     LS,
     *  ISRCHR,   RSRCHR,   IPT,      IELT,     IERR    )

C-----------------------------------------------------------------------

C DESCRIPTION:

C THIS SUBROUTINE CALCULATES THE CLOSEST POINT PROBLEM
C BETWEEN 'KOUNTS' PAIRS OF POINTS AND SURFACES.

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
C  P     LS       I/U    CONNECTIVITY OF ELEMENTS (4*NFSRF),
C                        NUMBERS REFER TO LOCATIONS IN XYZE ARRAY
C  P     ISRCHR   I/O    INTEGER SEARCH RESULTS
C  P     RSRCHR   I/O    REAL SEARCH RESULTS
C  P     IPT      I/U    POINT PAIRED WITH SURFACE LISTED IN IELT
C  P     IELT     I/U    SURFACE PAIRED WITH POINT LISTED IN IPT

C-----------------------------------------------------------------------

      include 'amesh.blk'
      include 'ebbyeb.blk'
      include 'toldat.blk'
      include 'tapes.blk'

C INPUT/OUTPUT ARRAYS
      DIMENSION
     *  XYZP(NPTS,NDIM)     ,XYZE(NPSRF,NDIM)  ,LS(NELNDA,NFSRF)    ,
     *  ISRCHR(NISR,NPTS)   ,RSRCHR(NRSR,NPTS)
      DIMENSION XX(27), YY(27), ZZ(27)

      IF( NISR .LT. 1 .OR. NRSR .LT. 3 .OR. NRSS .LT. 10 )THEN
        IERR = 1
        RETURN
      ENDIF

C check for Mesh-B point coincident with node of element in Mesh-A

      SIDE1 = (XYZE(LS(1,IELT),1)-XYZE(LS(2,IELT),1))**2
     &      + (XYZE(LS(1,IELT),2)-XYZE(LS(2,IELT),2))**2
      SIDE2 = (XYZE(LS(2,IELT),1)-XYZE(LS(3,IELT),1))**2
     &      + (XYZE(LS(2,IELT),2)-XYZE(LS(3,IELT),2))**2
      SIDE3 = (XYZE(LS(3,IELT),1)-XYZE(LS(4,IELT),1))**2
     &      + (XYZE(LS(3,IELT),2)-XYZE(LS(4,IELT),2))**2
      SIDE4 = (XYZE(LS(4,IELT),1)-XYZE(LS(1,IELT),1))**2
     &      + (XYZE(LS(4,IELT),2)-XYZE(LS(1,IELT),2))**2
      SIDMIN = MIN(SIDE1,SIDE2,SIDE3,SIDE4)
      SIDMAX = MAX(SIDE1,SIDE2,SIDE3,SIDE4)
      COTEST = EPS*EPS*SIDMIN
      DO 110 I = 1, 4
        A = XYZE(LS(I,IELT),1) - XYZP(IPT,1)
        B = XYZE(LS(I,IELT),2) - XYZP(IPT,2)
        DIST = A**2+B**2
        IF (DIST .LT. COTEST)THEN

C coincident node, so fill search results arrays
C no need to check for better search result

          INODE = I
          ISRCHR(1,IPT) = IELT
          CALL NODE (3,INODE,RSRCHR(1,IPT),RSRCHR(2,IPT),
     &      RSRCHR(3,IPT))
          GO TO 100
        END IF
 110  CONTINUE

C Mesh-B point not coincident with Mesh-A node so compute isoparametric
C coordinates. Use Newton's method

      SG = 0.
      TG = 0.
      RG = 0.
      ITER = 0

C Build Jacobian and invert

      DO 120 I = 1, NELNDA
        XX(I) = XYZE(LS(I,IELT),1)
        YY(I) = XYZE(LS(I,IELT),2)
        ZZ(I) = 0.
 120  CONTINUE
 130  CONTINUE
      CALL JACOBN (ITYPE,XX,YY,ZZ,SG,TG,RG,A11,A12,A13,A21,A22,A23,
     &  A31,A32,A33,F1,F2,F3)
      DETA = A11*A22 - A12*A21
      IF (ABS(DETA) .GT. 1.E-25)THEN

        AI11 =  A22/DETA
        AI12 = -A12/DETA
        AI21 = -A21/DETA
        AI22 =  A11/DETA

        FS = F1 - XYZP(IPT,1)
        FT = F2 - XYZP(IPT,2)
        SNEW = SG - (AI11*FS + AI12*FT)
        TNEW = TG - (AI21*FS + AI22*FT)

        ITER = ITER + 1
        DS = ABS(SNEW-SG)
        DT = ABS(TNEW-TG)
        IF (DS .LT. TOL .AND. DT .LT. TOL) GO TO 300
        SG = SNEW
        TG = TNEW
        IF (ITER .EQ. ITERMX)GO TO 100
        GO TO 130
      ELSE

C Zero Jacobian - check for degenerate quad (triangular element)

        TRITST = EPS*EPS*SIDMAX
        IF (SIDE1 .LT. TRITST)THEN
          XX(1) = XYZE(LS(1,IELT),1)
          XX(2) = XYZE(LS(3,IELT),1)
          XX(3) = XYZE(LS(4,IELT),1)
          YY(1) = XYZE(LS(1,IELT),2)
          YY(2) = XYZE(LS(3,IELT),2)
          YY(3) = XYZE(LS(4,IELT),2)
        ELSE IF (SIDE2 .LT. TRITST)THEN
          XX(1) = XYZE(LS(1,IELT),1)
          XX(2) = XYZE(LS(2,IELT),1)
          XX(3) = XYZE(LS(4,IELT),1)
          YY(1) = XYZE(LS(1,IELT),2)
          YY(2) = XYZE(LS(2,IELT),2)
          YY(3) = XYZE(LS(4,IELT),2)
        ELSE IF (SIDE3 .LT. TRITST)THEN
          XX(1) = XYZE(LS(1,IELT),1)
          XX(2) = XYZE(LS(2,IELT),1)
          XX(3) = XYZE(LS(3,IELT),1)
          YY(1) = XYZE(LS(1,IELT),2)
          YY(2) = XYZE(LS(2,IELT),2)
          YY(3) = XYZE(LS(3,IELT),2)
        ELSE IF (SIDE4 .LT. TRITST)THEN
          XX(1) = XYZE(LS(2,IELT),1)
          XX(2) = XYZE(LS(3,IELT),1)
          XX(3) = XYZE(LS(4,IELT),1)
          YY(1) = XYZE(LS(2,IELT),2)
          YY(2) = XYZE(LS(3,IELT),2)
          YY(3) = XYZE(LS(4,IELT),2)
        ELSE
          CALL ERROR ('QADSRC',
     &      'ZERO JACOBIAN FOUND DURING NEWTON ITERATION',
     &      'MESH-A ELEMENT',IELT,
     &      'ELEMENT IS NOT A DEGENERATE QUAD - GIVING UP',
     &      0,' ',' ',0)
          GO TO 100
        END IF

C Process as triangle

 210    CONTINUE
        CALL JACOBN (1,XX,YY,ZZ,SG,TG,RG,A11,A12,A13,A21,A22,A23,
     &    A31,A32,A33,F1,F2,F3)
        DETA = A11*A22 - A12*A21
        IF (ABS(DETA) .LT. 1.E-25)THEN
          CALL ERROR ('SRCHQ',
     &      'ZERO JACOBIAN FOUND DURING NEWTON ITERATION',
     &      'MESH-A ELEMENT',IELT,
     &      'TRYING TO PROCESS AS A DEGENERATE QUAD (TRIANGLE)',
     &      0,' ',' ',0)
        END IF

        AI11 =  A22/DETA
        AI12 = -A12/DETA
        AI21 = -A21/DETA
        AI22 =  A11/DETA

        FS = F1 - XYZP(IPT,1)
        FT = F2 - XYZP(IPT,2)
        SNEW = SG - (AI11*FS + AI12*FT)
        TNEW = TG - (AI21*FS + AI22*FT)

        ITER = ITER + 1
        DS = ABS(SNEW-SG)
        DT = ABS(TNEW-TG)
        IF (DS .LT. TOL .AND. DT .LT. TOL) GO TO 300
        SG = SNEW
        TG = TNEW
        IF (ITER .EQ. ITERMX)GO TO 100
        GO TO 210
      END IF

 300  CONTINUE

C Newton converged, load up search results arrays if appropriate

      IF (ABS(SNEW) .LT. STRLMT .AND. ABS(TNEW) .LT. STRLMT)THEN

C Search was adequate

        FTEST = MAX(ABS(RSRCHR(1,IPT)),ABS(RSRCHR(2,IPT)))
        FCOMP = MAX(ABS(SNEW),ABS(TNEW))
        IF (FTEST .GT. FCOMP .OR. ISRCHR(1,IPT) .EQ. 0)THEN

C New search is better, replace search results

          ISRCHR(1,IPT) = IELT
          RSRCHR(1,IPT) = SNEW
          RSRCHR(2,IPT) = TNEW
          RSRCHR(3,IPT) = 0.
        END IF
      END IF
 100  CONTINUE
      RETURN
      END
