C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      SUBROUTINE MKRNK(N,NTOTAL,NDIM,X,IND,IRNK,IRNK2)

C***********************************************************************

C     DESCRIPTION: THIS ROUTINE CONVERTS THE IRNK VECTORS IN THE SWEGLE
C                  SEARCH FROM AN INDIRECT GATHER TO A DIRECT GATHER

C       N        INTEGER   NUMBER OF ENTITIES THAT WAS SORTED
C       NDIM     INTEGER   NUMBER OF DIMENSIONS
C       X        REAL      ENTITIES TO BE SORTED

C       IND      INTEGER   INDEX VECTOR
C       IRNK     INTEGER   RANK VECTOR (INDIRECT)
C       IRNK2    INTEGER   RANK VECTOR (DIRECT)

C***********************************************************************

      include 'tapes.blk'

      DIMENSION X(NTOTAL,NDIM),IND(N,NDIM)
      DIMENSION IRNK(N,NDIM),IRNK2(N,NDIM,*)

       DO 11 IDM = 1, NDIM
         CALL INDEXX(X(1,IDM),IND(1,IDM),N,.true.)
         CALL RANK(N,IND(1,IDM),IRNK(1,IDM),N)
 11    CONTINUE

C CONSTRUCT DIRECT LISTS INTO ORDERED LIST OF POINTS
        IF(NDIM .EQ. 1)THEN
          DO 113 I = 1, N
            IRNK2(I,1,1) = IRNK(I,1)
113       CONTINUE
        ELSE IF( NDIM .EQ. 2)THEN
          DO 213 I = 1, N
            IRNK2(I,1,1) = IRNK(IND(I,1),2)
            IRNK2(I,2,1) = IRNK(IND(I,2),1)
213       CONTINUE
        ELSE IF( NDIM .EQ. 3)THEN
        DO 313 I=1,N
          IRNK2(I,1,1) = IRNK(IND(I,1),2)
          IRNK2(I,1,2) = IRNK(IND(I,1),3)
          IRNK2(I,2,1) = IRNK(IND(I,2),1)
          IRNK2(I,2,2) = IRNK(IND(I,2),3)
          IRNK2(I,3,1) = IRNK(IND(I,3),1)
          IRNK2(I,3,2) = IRNK(IND(I,3),2)
 313    CONTINUE
       ELSE
         PRINT*,'WRONG NUMBER OF DIMENSIONS IN MKRNK'
         STOP
       ENDIF

      RETURN
      END
