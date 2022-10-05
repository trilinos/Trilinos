C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      SUBROUTINE MKLSTV( NUMPTS,IND,IRNK2,IUP,ILO,INDX,
     *                   IE,LIST,NLIST,NBLKSZ,NSPC)

C-----------------------------------------------------------------------

C DESCRIPTION:

C VECTOR MAKE LIST (3D)
C GIVEN A LIST OF PARTICLES (IE  THEIR INDEX AND RANK) FIND
C THE LIST OF PARTICLES WITHIN THE BOUNDS SET BY XMIN AND XMAX
C FOR THE IE'TH PARTICLE IN THE VECTOR BLOCK

C-----------------------------------------------------------------------

C  CALLING ARGUMENTS:

C     NUMPTS INTEGER   NUMBER OF POINTS TO BE SEARCHED
C     IND    INTEGER   ORDER INDEX
C     IRNK2  INTEGER   RANK
C     IUP    INTEGER   SCRACTH (NBLKSZ LONG)
C     ILO    INTEGER   SCRACTH (NBLKSZ LONG)
C     INDX   INTEGER   SCRATCH (NBLKSZ LONG)
C     IE     INTEGER   PARTICLE NUMBER
C     LIST   INTEGER   LIST OF FOUND PARTICLES
C     NLIST  INTEGER   NUMBER OF PARTICLES FOUND
C     NBLKSZ INTEGER   BLOCK SIZE OF IUP AND ILO BLOCKS
C     NSPC   INTEGER   NUMBER OF SPACIAL COORD. (NUMBER OF DIMENSIONS)

C-----------------------------------------------------------------------

      DIMENSION
     *  IND(NUMPTS,NSPC),IRNK2(NUMPTS,NSPC,*),
     *  IUP(NBLKSZ,NSPC),INDX(NUMPTS),
     *  ILO(NBLKSZ,NSPC), LIST(NUMPTS)

C BUILD A LIST OF POINTS THAT ARE CLOSE TO SURFACE IE
      J = IE
      NLIST = 0
      IF( NSPC .EQ. 1)THEN
C============================== o n e   -  d ======================
        NUM1 = IUP(J,1) - ILO(J,1) + 1
        ILOW = ILO(J,1)
        IUPR = IUP(J,1)
        DO 101 I1 = ILOW, IUPR
          NLIST = NLIST +1
          LIST(NLIST) = IND(I1,1)
 101    CONTINUE

      ELSE IF( NSPC .EQ. 2 )THEN
C============================== t w o   -  d ======================
        NUM1 = IUP(J,1) - ILO(J,1) + 1
        NUM2 = IUP(J,2) - ILO(J,2) + 1
C DO WE HAVE A LIST ?
        IF( NUM2.LE.0 .OR. NUM1.LE.0 ) RETURN
C SELECT WHICH LIST IS THE SMALLEST
        IF( NUM1 .LE. NUM2 )THEN
          IXYZ = 1
          IY = 2
          NUM = NUM1
        ELSE
          IXYZ = 2
          IY = 1
          NUM = NUM2
        ENDIF

        ILOW = ILO(J,IXYZ)
        IUPR = IUP(J,IXYZ)
C FIRST TEST
        IF( NUM .GT. 64 ) THEN
          DO 201 I1 = ILOW, IUPR
            IF( IRNK2(I1,IXYZ,1) .GE. ILO(J,IY) .AND.
     *        IRNK2(I1,IXYZ,1) .LE. IUP(J,IY) )THEN
              NLIST = NLIST +1
              LIST(NLIST) = IND(I1,IXYZ)
            ENDIF
 201      CONTINUE
        ELSE
          DO 202 I1 = ILOW, IUPR
            IF( IRNK2(I1,IXYZ,1) .GE. ILO(J,IY) .AND.
     *        IRNK2(I1,IXYZ,1) .LE. IUP(J,IY) )THEN
              NLIST = NLIST +1
              LIST(NLIST) = IND(I1,IXYZ)
            ENDIF
 202      CONTINUE
        ENDIF

      ELSE IF( NSPC .EQ. 3 )THEN
C============================== t h r e e   -  d ======================
        NUM1 = IUP(J,1) - ILO(J,1) + 1
        NUM2 = IUP(J,2) - ILO(J,2) + 1
        NUM3 = IUP(J,3) - ILO(J,3) + 1
C DO WE HAVE A LIST ?
        IF( NUM3.LE.0 .OR. NUM2.LE.0 .OR. NUM1.LE.0 ) RETURN
C SELECT WHICH LIST IS THE SMALLEST
        IF( NUM1 .LE. NUM2 .AND. NUM1 .LE. NUM3 )THEN
          IXYZ = 1
          IY = 2
          IZ = 3
          NUM = NUM1
        ELSEIF( NUM2 .LE. NUM1 .AND. NUM2 .LE. NUM3 )THEN
          IXYZ = 2
          IY = 1
          IZ = 3
          NUM = NUM2
        ELSE
          IXYZ = 3
          IY = 1
          IZ = 2
          NUM = NUM3
        ENDIF

        ILOW = ILO(J,IXYZ)
        IUPR = IUP(J,IXYZ)
        IF (ILOW.EQ.0) THEN
          NLIST = 0
          RETURN
        ENDIF
        ILP = 0
C FIRST TEST
        IF( NUM .GT. 64 ) THEN
          DO 301 I1 = ILOW, IUPR
            IF( IRNK2(I1,IXYZ,1) .GE. ILO(J,IY) .AND.
     *        IRNK2(I1,IXYZ,1) .LE. IUP(J,IY) )THEN
              ILP = ILP +1
              INDX(ILP) = I1
            ENDIF
 301      CONTINUE
        ELSE
          DO 302 I1 = ILOW, IUPR
            IF( IRNK2(I1,IXYZ,1) .GE. ILO(J,IY) .AND.
     *        IRNK2(I1,IXYZ,1) .LE. IUP(J,IY) )THEN
              ILP = ILP +1
              INDX(ILP) = I1
            ENDIF
 302      CONTINUE
        ENDIF
C SECOND TEST
        IF( ILP .GT. 64 ) THEN
          DO 311 I1 = 1, ILP
            IF( IRNK2(INDX(I1),IXYZ,2) .GE. ILO(J,IZ) .AND.
     *        IRNK2(INDX(I1),IXYZ,2) .LE. IUP(J,IZ) )THEN
              NLIST = NLIST + 1
              LIST(NLIST) = IND(INDX(I1),IXYZ)
            ENDIF
 311      CONTINUE
        ELSE
          DO 313 I1 = 1, ILP
            IF( IRNK2(INDX(I1),IXYZ,2) .GE. ILO(J,IZ) .AND.
     *        IRNK2(INDX(I1),IXYZ,2) .LE. IUP(J,IZ) )THEN
              NLIST = NLIST + 1
              LIST(NLIST) = IND(INDX(I1),IXYZ)
            ENDIF
 313      CONTINUE
        ENDIF
      ENDIF

      RETURN
      END

