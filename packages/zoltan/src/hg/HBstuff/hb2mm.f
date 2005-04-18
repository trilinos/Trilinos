c     *********** Harwell-Boeing to Matrix Market format ***********
c
c     This program reads a sparse matrix in Harwell-Boeing format from input
c     and writes it in Matrix Market format without headers to output. 
c     The output indices start numbering from 1.
c     See Harwell-Boeing users guide.
c
c     Original version by Ann Sweeney.
c     Updated by Rob Bisseling and Erik Boman.

      INTEGER    MAXSIZ
c     MAXSIZ is maximum no. of nonzeros - Adjust as necessary!
      PARAMETER (MAXSIZ=25000000)
      PARAMETER (MAXrhs=100000)
      INTEGER    LUNIT,    LOUT
      PARAMETER (LUNIT=5, LOUT=6)
      INTEGER I, J, COLIND(MAXSIZ)

      CHARACTER TITLE*72, KEY*8, MXTYPE*3, RHSTYP*3,
     :          PTRFMT*16, INDFMT*16, VALFMT*20, RHSFMT*20
      INTEGER TOTCRD, PTRCRD, INDCRD, VALCRD, RHSCRD, 
     :        NROW, NCOL, NNZERO, NELTVL, NRHS, NRHSIX, NRHSVL,
     :        NGUESS, NEXACT
      INTEGER POINTR(MAXSIZ), ROWIND(MAXSIZ), RHSPTR(MAXrhs), 
     :        RHSIND(MAXrhs)
      REAL VALUES(MAXSIZ), RHSVAL(MAXrhs), XEXACT(MAXrhs),
     :        SGUESS(MAXrhs)

c    *************** read in header block ***************
      READ (LUNIT,1000) TITLE, KEY, TOTCRD, PTRCRD, INDCRD, VALCRD, 
     :                  RHSCRD, MXTYPE, NROW, NCOL, NNZERO, NELTVL,
     :                  PTRFMT, INDFMT, VALFMT, RHSFMT

      IF (RHSCRD.GT.0) READ (LUNIT, 1001) RHSTYP, NRHS, NRHSIX

 1000 FORMAT (A72,A8/5I14/A3,11X,4I14/2A16,2A20)
 1001 FORMAT (A3,11X,2I14)
     
c    *************** read in matrix structure ***************
      READ (LUNIT, PTRFMT) (POINTR(I),I=1,NCOL+1)
      READ (LUNIT, INDFMT) (ROWIND(I),I=1,NNZERO) 
      IF (VALCRD.GT.0) THEN
c         ******* read matrix values *****
          IF (MXTYPE(3:3).EQ.'A') THEN
              READ(LUNIT,VALFMT) (VALUES(I), I=1,NNZERO)
          ELSE
              READ(LUNIT,VALFMT) (VALUES(I), I=1,NELTVL)
          ENDIF
c         ******* read right-hand sides *****
          IF (NRHS.GT.0) THEN
              IF (RHSTYP(1:1).EQ.'F') THEN
c                 ******* read dense right-hand sides *******
                  NRHSVL = NROW*NRHS
                  READ(LUNIT,RHSFMT) (RHSVAL(I),I=1,NRHSVL)
              ELSE
c                 ******* read sparse or elemental right-hand sides 
                  IF (MXTYPE(3:3).EQ.'A') THEN
c                     *** sparse right-hand sides: read pointer array
                      READ (LUNIT,PTRFMT) (RHSPTR (I),I=1,NRHS+1)
c                     *** read sparsity pattern for right-hand sides
                      READ (LUNIT,INDFMT) (RHSIND(I),I=1,NRHSIX)
c                     *** read sparse right-hand values 
                      READ (LUNIT,RHSFMT) (RHSVAL(I),I=1,NRHSIX)
                  ELSE
c                     *** read elemental right-hand sides
                      NRHSVL = NNZERO * NRHS
                      READ (LUNIT,RHSFMT) (RHSVAL(I),I=1,NRHSVL)
                  ENDIF
              ENDIF
              IF (RHSTYP(2:2).EQ.'G') THEN
c                 ******* read starting guesses
                  NGUESS = NROW * NRHS
                  READ (LUNIT,RHSFMT) (SGUESS(I),I=1,NGUESS)
              ENDIF
              IF (RHSTYP(3:3).EQ.'X') THEN
c                 ******* read solution vectors
                  NEXACT = NROW * NRHS
                  READ (LUNIT,RHSFMT) (XEXACT(I),I=1,NEXACT)
              ENDIF
          ENDIF
      ENDIF
 
c     **********  write on output file in matrix market format ************
c     **********  numbered from 1                              ************

      WRITE(LOUT,*) '%%MatrixMarket matrix coordinate real general'
      WRITE(LOUT,*) '% ', TITLE
      WRITE(LOUT,*) NROW,NCOL,NNZERO
      DO 20 J=1,NCOL
          DO 10 I=POINTR(J),POINTR(J+1)-1
              COLIND(I) = J
   10     CONTINUE
   20 CONTINUE

      DO 25 J=1,NNZERO
          WRITE(LOUT, *) ROWIND(J), COLIND(J), VALUES(J)
   25 CONTINUE
  
      STOP
      END
