C-------------------------------------------------------------------------------
C 1.  toms/575
C keywords:  nonsymmetric permutations, maximum transversal, maximum assignment
C gams:  D2e
C title:  MC21A
C for:  row permutation for a zero-free diagonal That is, given the pattern of
C nonzeros of a sparse matrix, this routine attempts to find a permutation of
C its rows that makes the matrix have no zeros on its diagonal
C by:  I.S. Duff
C ref:  ACM TOMS 7 (1981) 387-390
C-------------------------------------------------------------------------------

C S IS THE STANDARD FORTRAN VERSION   SI/                                      1
C I IS THE IBM VERSION                                                         2
C                                                                              3
      SUBROUTINE MC21A(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW)
C
C DESCRIPTION OF PARAMETERS.
C INPUT VARIABLES   N,ICN,LICN,IP,LENR.
C OUTPUT VARIABLES  IPERM,NUMNZ.
C
C N   ORDER OF MATRIX.
C ICN ARRAY CONTAINING THE COLUMN INDICES OF THE NON-ZEROS.  THOSE
C     BELONGING TO A SINGLE ROW MUST BE CONTIGUOUS BUT THE ORDERING
C     OF COLUMN INDICES WITHIN EACH ROW IS UNIMPORTANT AND WASTED
C     SPACE BETWEEN ROWS IS PERMITTED.
C LICN  LENGTH OF ARRAY ICN.
C IP  IP(I), I=1,2,...N, IS THE POSITION IN ARRAY ICN OF THE FIRST
C     COLUMN INDEX OF A NON-ZERO IN ROW I.
C LENR  LENR(I) IS THE NUMBER OF NON-ZEROS IN ROW I, I=1,2,..N.
C IPERM CONTAINS PERMUTATION TO MAKE DIAGONAL HAVE THE SMALLEST
C     NUMBER OF ZEROS ON IT.  ELEMENTS (IPERM(I),I) I=1, ... N ARE
C     NON-ZERO AT THE END OF THE ALGORITHM UNLESS MATRIX
C     IS STRUCTURALLY SINGULAR.  IN THIS CASE, (IPERM(I),I) WILL
C     BE ZERO FOR N-NUMNZ ENTRIES.
C NUMNZ NUMBER OF NON-ZEROS ON DIAGONAL OF PERMUTED MATRIX.
C IW  WORK ARRAY  ..  SEE LATER COMMENTS.
C
      INTEGER IP(N)
C     INTEGER*2 ICN(LICN),LENR(N),IPERM(N),IW(N,4)    I/
      INTEGER   ICN(LICN),LENR(N),IPERM(N),IW(N,4)
      CALL MC21B(N,ICN,LICN,IP,LENR,IPERM,NUMNZ,IW(1,1),IW(1,2),IW(1,3),
     1IW(1,4))
      RETURN
      END

