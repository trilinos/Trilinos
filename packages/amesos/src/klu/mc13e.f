C-------------------------------------------------------------------------------
C 2.  toms/529  
C keywords: symmetric permutations, block triangular, depth first search, sparse
C gams:  D2e  
C title:  MC13D  for:  finding symmetric permutations to block triangular form.
C That is, given the column numbers of the nonzeros in each row of a sparse
C matrix, this subroutine finds a symmetric permutation that makes the matrix
C block lower triangular.
C by:  I.S. Duff and J.K. Reid
C ref:  ACM TOMS 4 (1978) 189-192
C-------------------------------------------------------------------------------

      SUBROUTINE MC13E(N, ICN, LICN, IP, LENR, ARP, IB, NUM, LOWL,
     *  NUMB, PREV)
      INTEGER STP, DUMMY
      INTEGER IP(N)
C
C ARP(I) IS ONE LESS THAN THE NUMBER OF UNSEARCHED EDGES LEAVING
C     NODE I.  AT THE END OF THE ALGORITHM IT IS SET TO A
C     PERMUTATION WHICH PUTS THE MATRIX IN BLOCK LOWER
C     TRIANGULAR FORM.
C IB(I) IS THE POSITION IN THE ORDERING OF THE START OF THE ITH
C     BLOCK.  IB(N+1-I) HOLDS THE NODE NUMBER OF THE ITH NODE
C     ON THE STACK.
C LOWL(I) IS THE SMALLEST STACK POSITION OF ANY NODE TO WHICH A PATH
C     FROM NODE I HAS BEEN FOUND.  IT IS SET TO N+1 WHEN NODE I
C     IS REMOVED FROM THE STACK.
C NUMB(I) IS THE POSITION OF NODE I IN THE STACK IF IT IS ON
C     IT, IS THE PERMUTED ORDER OF NODE I FOR THOSE NODES
C     WHOSE FINAL POSITION HAS BEEN FOUND AND IS OTHERWISE ZERO.
C PREV(I) IS THE NODE AT THE END OF THE PATH WHEN NODE I WAS
C     PLACED ON THE STACK.
C     INTEGER*2 ICN(LICN),LENR(N),ARP(N),IB(N),LOWL(N),NUMB(N),       I/
C    1PREV(N)                                                         I/
      INTEGER ICN(LICN), LENR(N), ARP(N), IB(N), LOWL(N), NUMB(N),
     *  PREV(N)
C
C
C   ICNT IS THE NUMBER OF NODES WHOSE POSITIONS IN FINAL ORDERING HAVE
C     BEEN FOUND.
      ICNT = 0
C NUM IS THE NUMBER OF BLOCKS THAT HAVE BEEN FOUND.
      NUM = 0
      NNM1 = N + N - 1
C
C INITIALIZATION OF ARRAYS.
      DO 10 J=1,N
        NUMB(J) = 0
        ARP(J) = LENR(J) - 1
   10 CONTINUE
C
C
      DO 90 ISN=1,N
C LOOK FOR A STARTING NODE
        IF (NUMB(ISN).NE.0) GO TO 90
        IV = ISN
C IST IS THE NUMBER OF NODES ON THE STACK ... IT IS THE STACK POINTER.
        IST = 1
C PUT NODE IV AT BEGINNING OF STACK.
        LOWL(IV) = 1
        NUMB(IV) = 1
        IB(N) = IV
C
C THE BODY OF THIS LOOP PUTS A NEW NODE ON THE STACK OR BACKTRACKS.
        DO 80 DUMMY=1,NNM1
          I1 = ARP(IV)
C HAVE ALL EDGES LEAVING NODE IV BEEN SEARCHED.
          IF (I1.LT.0) GO TO 30
          I2 = IP(IV) + LENR(IV) - 1
          I1 = I2 - I1
C
C LOOK AT EDGES LEAVING NODE IV UNTIL ONE ENTERS A NEW NODE OR
C     ALL EDGES ARE EXHAUSTED.
          DO 20 II=I1,I2
            IW = ICN(II)
C HAS NODE IW BEEN ON STACK ALREADY.
            IF (NUMB(IW).EQ.0) GO TO 70
C UPDATE VALUE OF LOWL(IV) IF NECESSARY.
            IF (LOWL(IW).LT.LOWL(IV)) LOWL(IV) = LOWL(IW)
   20     CONTINUE
C
C THERE ARE NO MORE EDGES LEAVING NODE IV.
          ARP(IV) = -1
C IS NODE IV THE ROOT OF A BLOCK.
   30     IF (LOWL(IV).LT.NUMB(IV)) GO TO 60
C
C ORDER NODES IN A BLOCK.
          NUM = NUM + 1
          IST1 = N + 1 - IST
          LCNT = ICNT + 1
C PEEL BLOCK OFF THE TOP OF THE STACK STARTING AT THE TOP AND
C     WORKING DOWN TO THE ROOT OF THE BLOCK.
          DO 40 STP=IST1,N
            IW = IB(STP)
            LOWL(IW) = N + 1
            ICNT = ICNT + 1
            NUMB(IW) = ICNT
            IF (IW.EQ.IV) GO TO 50
   40     CONTINUE
   50     IST = N - STP
          IB(NUM) = LCNT
C ARE THERE ANY NODES LEFT ON THE STACK.
          IF (IST.NE.0) GO TO 60
C HAVE ALL THE NODES BEEN ORDERED.
          IF (ICNT.LT.N) GO TO 90
          GO TO 100
C
C BACKTRACK TO PREVIOUS NODE ON PATH.
   60     IW = IV
          IV = PREV(IV)
C UPDATE VALUE OF LOWL(IV) IF NECESSARY.
          IF (LOWL(IW).LT.LOWL(IV)) LOWL(IV) = LOWL(IW)
          GO TO 80
C
C PUT NEW NODE ON THE STACK.
   70     ARP(IV) = I2 - II - 1
          PREV(IW) = IV
          IV = IW
          IST = IST + 1
          LOWL(IV) = IST
          NUMB(IV) = IST
          K = N + 1 - IST
          IB(K) = IV
   80   CONTINUE
C
   90 CONTINUE
C
C
C PUT PERMUTATION IN THE REQUIRED FORM.
  100 DO 110 I=1,N
        II = NUMB(I)
        ARP(II) = I
  110 CONTINUE
      RETURN
      END

