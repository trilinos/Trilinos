c call MC21B (n, II, nz, Cp, W, Zperm, ndiag, Offp, Cperm, Pr, Pc)
c	   MC21B calling interface:
c	   input:	n, II (1..nz), nz, Cp (n), W (n):
c			n-by-n matrix, col is of length W (col),
c			and its pattern is located in
c			II (Cp (col) ... Cp (col)+W(col)-1)
c	   output:	Zperm (n), the permutation, such that
c			colold = Zperm (col), and ndiag (number of 
c			structural nonzeros on the diagonal.	
c	   		matrix is structurally singular if ndiag < n
c	   workspace:	Offp, Cperm, Pr, Pc


C/                                                                            34
      SUBROUTINE MC21B(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,PR,ARP,CV,OUT)
      INTEGER IP(N)
C
C     DIVISION OF WORK ARRAY IS NOW DESCRIBED.
C
C PR(I) IS THE PREVIOUS ROW TO I IN THE DEPTH FIRST SEARCH.
C ARP(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I
C     WHICH HAVE NOT BEEN SCANNED WHEN LOOKING FOR A CHEAP ASSIGNMENT.
C CV(I) IS THE MOST RECENT ROW EXTENSION AT WHICH COLUMN I
C     WAS VISITED.
C OUT(I) IS ONE LESS THAN THE NUMBER OF NON-ZEROS IN ROW I
C     WHICH HAVE NOT BEEN SCANNED DURING ONE PASS THROUGH THE
C     MAIN LOOP.
C
C     INTEGER*2 ICN(LICN),LENR(N),IPERM(N),PR(N),CV(N),  I/
C    1ARP(N),OUT(N)  I/
      INTEGER   ICN(LICN),LENR(N),IPERM(N),PR(N),CV(N),
     1ARP(N),OUT(N)
C
C   INITIALIZATION OF ARRAYS.
      DO 10 I=1,N
      ARP(I)=LENR(I)-1
      CV(I)=0
      IPERM(I)=0
   10 continue
      NUMNZ=0
C
C
C   MAIN LOOP.
C   EACH PASS ROUND THIS LOOP EITHER RESULTS IN A NEW ASSIGNMENT
C OR GIVES A ROW WITH NO ASSIGNMENT.
      DO 130 JORD=1,N
      J=JORD
      PR(J)=-1
      DO 100 K=1,JORD
C LOOK FOR A CHEAP ASSIGNMENT
      IN1=ARP(J)
      IF (IN1.LT.0) GO TO 60
      IN2=IP(J)+LENR(J)-1
      IN1=IN2-IN1
      DO 50 II=IN1,IN2
      I=ICN(II)
      IF (IPERM(I).EQ.0) GO TO 110
   50 CONTINUE
C   NO CHEAP ASSIGNMENT IN ROW.
      ARP(J)=-1
C   BEGIN LOOKING FOR ASSIGNMENT CHAIN STARTING WITH ROW J.
   60 OUT(J)=LENR(J)-1
C INNER LOOP.  EXTENDS CHAIN BY ONE OR BACKTRACKS.
      DO 90 KK=1,JORD
      IN1=OUT(J)
      IF (IN1.LT.0) GO TO 80
      IN2=IP(J)+LENR(J)-1
      IN1=IN2-IN1
C FORWARD SCAN.
      DO 70 II=IN1,IN2
      I=ICN(II)
      IF (CV(I).EQ.JORD) GO TO 70
C   COLUMN I HAS NOT YET BEEN ACCESSED DURING THIS PASS.
      J1=J
      J=IPERM(I)
      CV(I)=JORD
      PR(J)=J1
      OUT(J1)=IN2-II-1
      GO TO 100
   70 CONTINUE
C
C   BACKTRACKING STEP.
   80 J=PR(J)
      IF (J.EQ.-1) GO TO 130
   90 CONTINUE
C
  100 CONTINUE
C
C   NEW ASSIGNMENT IS MADE.
  110 IPERM(I)=J

      ARP(J)=IN2-II-1
      NUMNZ=NUMNZ+1
      DO 120 K=1,JORD
      J=PR(J)
      IF (J.EQ.-1) GO TO 130
      II=IP(J)+LENR(J)-OUT(J)-2
      I=ICN(II)
      IPERM(I)=J
  120 CONTINUE
C
  130 CONTINUE
C
C   IF MATRIX IS STRUCTURALLY SINGULAR, WE NOW COMPLETE THE
C PERMUTATION IPERM.
      IF (NUMNZ.EQ.N) GO TO 500
      DO 140 I=1,N
      ARP(I)=0
  140 continue
      K=0
      DO 160 I=1,N
      IF (IPERM(I).NE.0) GO TO 150
      K=K+1
      OUT(K)=I
      GO TO 160
  150 J=IPERM(I)
      ARP(J)=I
  160 CONTINUE
      K=0
      DO 170 I=1,N
      IF (ARP(I).NE.0) GO TO 170
      K=K+1
      IOUTK=OUT(K)
      IPERM(IOUTK)=I
  170 CONTINUE
  500 RETURN
      END

