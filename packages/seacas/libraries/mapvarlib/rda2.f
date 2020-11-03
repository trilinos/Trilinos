C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
*DECK,RDA2
      SUBROUTINE RDA2 (IDBLKA,ICONA,NDLSTA,STATUS,MAXLN)

C     ******************************************************************

C     SUBROUTINE TO READ MESH A, WRITE MESH C DATA AS APPROPRIATE

C     Calls subroutine ERROR

C     Called by MAPVAR

C     ******************************************************************

C  IDBLKA  INT  Element block I.D. donor mesh
C  ICONA   INT  Connectivity for elt block (1:nelnda,1:numeba)
C  NDLSTA  INT  The array that identifies the local element block node
C               number with the global mesh node number (1:numnda)
C  STATUS  REAL Element status - alive or dead (1:numeba)
C  ITYPE   INT  Element type
C  NELNDA  INT  Number of nodes per element
C  NUMNDA  INT  Number of nodes in element block
C  NUMEBA  INT  Number of elements in element block
C  MAXLN   INT  Maximum number of elements per node for INVCON

C     ******************************************************************

      CHARACTER*(32) TYP

      include 'aexds1.blk'
      include 'aexds2.blk'
      include 'amesh.blk'
      include 'ebbyeb.blk'
      include 'ex2tp.blk'
      include 'steps.blk'

      DIMENSION ICONA(NELNDA,*),NDLSTA(*),STATUS(*)

C     ******************************************************************

C element type per element block

C fix this routine when i have time
C create array 1-nnodes
C loop over all elements - add 1 to value in array whenever
C                          node appears in connectivity
C maxln = max value of array

        CALL EXGELB (NTP2EX,IDBLKA,TYP,NUMEBA,NELNDA,
     &               NATRIB,IERR)
        CALL EXUPCS(TYP)

        IF (TYP(1:3) .EQ. 'QUA')THEN
          IF (NELNDA .EQ. 4)THEN
            ITYPE = 3
          ELSE IF (NELNDA .EQ. 8)THEN
            ITYPE = 4
          ELSE IF (NELNDA .EQ. 9)THEN
            ITYPE = 5
          ELSE
            CALL ERROR ('RDA2','UNSUPPORTED ELEMENT TYPE',' ',0,' ',0,
     1                  'TYPE',typ,1)
          END IF
        ELSE IF (TYP(1:3) .EQ. 'HEX') THEN
          ITYPE = 10
        ELSE IF (TYP(1:3) .EQ. 'SHE')THEN
          ITYPE = 13
        ELSE IF (TYP(1:3) .EQ. 'TET') THEN
          ITYPE = 6
        ELSE
          CALL ERROR ('RDA2','UNSUPPORTED ELEMENT TYPE',' ',0,' ',0,
     1    'TYPE',typ,1)
        END IF

      CALL EXGELC(NTP2EX,IDBLKA,ICONA(1,1),IERR)

      DO 5 I = 1, NODESA
        NDLSTA(I) = 0
 5    CONTINUE

      DO 10 IEL = 1, NUMEBA
        DO 20 INODE = 1, NELNDA
          NDLSTA(ICONA(INODE,IEL)) = NDLSTA(ICONA(INODE,IEL)) + 1
 20     CONTINUE
 10   CONTINUE

      NUMNDA = 0

      MAXLN = 0
      DO 30 I = 1, NODESA
        IF (NDLSTA(I) .GT. 0) THEN
           if (ndlsta(i) .gt. maxln) then
              maxln = ndlsta(i)
           end if
          NUMNDA = NUMNDA + 1
          NDLSTA(NUMNDA) = I
        END IF
 30   CONTINUE

C get STATUS array for use in SEARCH so that dead elements can be
C eliminated from the search

      DO 99 I = 1, NUMEBA
         STATUS(I) = 0.
   99 CONTINUE
      DO 100 ISTATUS = 1, NVAREL
         IF (NAMVAR(nvargp+ISTATUS) .NE. 'STATUS')GO TO 100
         CALL EXGEV(NTP2EX,ISTEP,ISTATUS,IDBLKA,NUMEBA,STATUS,IERR)
  100 CONTINUE

      RETURN
      END
