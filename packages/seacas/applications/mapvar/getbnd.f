C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      SUBROUTINE GETBND(LBLK,NE,X,IND,NP,XMIN,XMAX,NDIM,ILO,IUP,
     *                  ISCR,RSCR )

C-----------------------------------------------------------------------

C DESCRIPTION:

C  FIND THE ELEMENTS IN A SORTED ARRAY X WHOSE VALUES FALL IN THE
C  INTERVAL BETWEEN XMIN AND XMAX. NO ELEMENTS HAVING
C  VALUES EQUAL TO XMIN OR XMAX ARE INCLUDED. SINCE THE ARRAY IS
C  SORTED, THE ELEMENTS CAN BE SPECIFIED BY THE UPPER AND
C  LOWER ELEMENT NUMBERS IN THE RANGE.

C     |    X(ILO)   . . . . . . . . . . . . .   X(IUP)    |
C    XMIN                                               XMAX       X>

C  IT IS ASSUMED THAT THE ARRAY X HAS BEEN SORTED IN INCREASING ORDER,
C  BUT THE ELEMENTS HAVE NOT BEEN MOVED.
C  THE SORTED LIST IS DETERMINED BY THE ARRAY INDX,
C  WHICH POSITIONS THE ORIGINAL UNSORTED X ARRAY ELEMENTS
C  IN THE SORTED LIST. THUS, THE 5TH ELEMENT IN THE SORTED LIST IS
C    X(IND(5))

C-----------------------------------------------------------------------

C  INPUT:

C  X      -  array in unsorted order
C  IND    -  index array giving the element order in the sorted list
C  NP     -  the number of particles in the list
C  XMIN   -  the lower limit of the interval
C  XMAX   -  the upper limit of the interval
C  NDIM   -  the dimension of the arrays

C  OUTPUT:

C  ILO    -  the first element in the sorted list .gt. xmin
C  IUP    -  the last element in the sorted list .lt. xmax

C-----------------------------------------------------------------------

      include 'tapes.blk'

      DIMENSION
     *  X(NDIM),IND(NDIM),XMIN(LBLK),XMAX(LBLK),ILO(LBLK),IUP(LBLK)
      DIMENSION
     *  ISCR(5*LBLK),RSCR(LBLK)

C INTEGER SCRATCH SPACE
      ISP    = 1
      LIL    = ISP
      ISP    = ISP + LBLK
      LIU    = ISP
      ISP    = ISP + LBLK
      LIT    = ISP
      ISP    = ISP + LBLK
      LINDX1 = ISP
C      ISP    = ISP + LBLK
C      LINDX2 = ISP
C      ISP    = ISP + LBLK
C REAL SCRATCH SPACE
      ISP   = 1
      LXTST = ISP
C      ISP   = ISP + LBLK

      DO 200 J = 1, NE, LBLK
            N = MIN(LBLK,NE-J+1)
C  SEARCH TO FIND THE FIRST ELEMENT .GE. XMIN
      CALL SRCHGE(LBLK,N,X,IND,XMIN(J),1,     NP,NDIM,ILO(J),
     *     ISCR(LIL),ISCR(LIU),ISCR(LIT),ISCR(LINDX1),ISCR(LINDX1),
     *     RSCR(LXTST) )
C  SEARCH TO FIND THE FIRST ELEMENT .GT. XMAX
      CALL SRCHGT(LBLK,N,X,IND,XMAX(J),ILO(J),NP,NDIM,IUP(J),
     *     ISCR(LIL),ISCR(LIU),ISCR(LIT),ISCR(LINDX1),ISCR(LINDX1),
     *     RSCR(LXTST) )
C  THE PREVIOUS ELEMENT IS THE LAST ONE .LT. XMAX
      DO 100 JJ = 1, N
       IUP(J+JJ-1)=IUP(J+JJ-1) - 1
 100  CONTINUE

 200  CONTINUE
      RETURN
      END

