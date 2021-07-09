C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      SUBROUTINE SRCHGE( LBLK,NE,X,IND,XV,IMIN,IMAX,NDIM,I,
     *                  IL,IU,IT,INDX1,INDX2,XTST )

C-----------------------------------------------------------------------

C DESCRIPTION:

C  PERFORM A BINARY SEARCH TO FIND THE ELEMENT NUMBER I
C  OF A SORTED ARRAY FOR WHICH ALL ELEMENTS AT I OR ABOVE ARE
C  GREATER OR EQUAL TO THAN SOME VALUE XV,
C  WHILE ALL ELEMENTS BELOW I ARE LESS THAN THE VALUE.

C       X(I-2)     X(I-1)    X(I)    X(I+1)   X(I+2)
C                      XV                             X>

C  ASSUMED THAT THE ARRAY X HAS BEEN SORTED IN INCREASING ORDER,
C  BUT THE ELEMENTS HAVE NOT BEEN MOVED.
C  THE SORTED LIST IS DETERMINED BY THE ARRAY INDX,
C  WHICH POSITIONS THE ORIGINAL UNSORTED X ARRAY ELEMENTS
C  IN THE SORTED LIST. THUS, THE 5TH ELEMENT IN THE SORTED LIST IS
C    X(IND(5))

C-----------------------------------------------------------------------

C  INPUT:

C  X      -  ARRAY IN UNSORTED ORDER
C  IND    -  INDEX ARRAY GIVING THE ELEMENT ORDER IN THE SORTED LIST
C  XV     -  X VALUE TO TEST AGAINST
C  IMIN   -  THE LOWEST NUMBERED POSITION IN THE SORTED LIST TO TEST
C  IMAX   -  THE HIGHEST NUMBERED POSITION IN THE SORTED LIST TO TEST
C  NDIM   -  THE DIMENSION OF THE ARRAYS

C  OUTPUT:

C  I      -  THE FIRST POSITION IN THE SORTED LIST .GT. XV

C  SCRATCH:

C  IL
C  IU
C  IT
C  INDX1
C  INDX2
C  XTST

C-----------------------------------------------------------------------

      DIMENSION
     *  X(NDIM), IND(NDIM), XV(LBLK), I(LBLK)
      DIMENSION
     *  IL(LBLK),IU(LBLK),IT(LBLK),INDX1(LBLK),INDX2(LBLK),
     *  XTST(LBLK)

      IF (IMAX.EQ.0.OR.NDIM.EQ.0) THEN
         DO J = 1, NE
            I(J) = 0
         ENDDO
         RETURN
      ENDIF
      DO 25 J = 1, NE
        IL(J) = IMIN
        IU(J) = IMAX
        INDX1(J) = J
 25   CONTINUE
        ILOOP = NE
 1000 CONTINUE

      DO 50 JJ = 1, ILOOP
        INDX2(JJ) = INDX1(JJ)
        J = INDX1(JJ)
        IT(J) =  (IU(J) + IL(J)) / 2
 50   CONTINUE
      DO 35 J = 1, NE
       XTST(J) = X( IND(IT(J)) )
 35   CONTINUE

      IF ( ILOOP .GT. 64) THEN

      ILP = 0
      DO 60 JJ = 1, ILOOP
        J = INDX2(JJ)
        IF( XTST(J) .LT. XV(J) )THEN
          IL(J) =IT(J) + 1
        ELSE
          IU(J) =IT(J) - 1
        ENDIF
        IF( IL(J) .LE. IU(J)) THEN
          ILP = ILP + 1
          INDX1(ILP) = J
        ENDIF
 60   CONTINUE
      ELSE

      ILP = 0
      DO 51 JJ = 1, ILOOP
        J = INDX2(JJ)
        IF( XTST(J) .LT. XV(J) )THEN
          IL(J) =IT(J) + 1
        ELSE
          IU(J) =IT(J) - 1
        ENDIF
        IF( IL(J) .LE. IU(J)) THEN
          ILP = ILP + 1
          INDX1(ILP) = J
        ENDIF
 51   CONTINUE
      ENDIF

      ILOOP = ILP

      IF(ILOOP .NE. 0 )GO TO 1000
C  RANGE HAD NARROWED TO 1 LOCATION. HOWEVER, THE POINT LAST TESTED
C  COULD BE ABOVE, BELOW, OR ON THE SEARCH POINT. CHECK FOR PROPER CASE
      DO 200 J = 1, NE
        IF( XTST(J) .LT. XV(J) )THEN
          I(J) = IT(J) + 1
        ELSE
          I(J) = IT(J)
        ENDIF
 200  CONTINUE

      RETURN
      END

