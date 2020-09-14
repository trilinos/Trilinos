C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
*DECK,RDB2
      SUBROUTINE RDB2 (IDBLKB,IDBLKA,ICONB,NDLSTB)

C     ******************************************************************

C     SUBROUTINE TO GET MESH B DATA

C     Calls subroutine ERROR

C     Called by MAPVAR

C     ******************************************************************

C  IDBLKA  INT   Element block I.D. - donor mesh
C  IDBLKB  INT   Element block I.D. - recipient mesh
C  ICONB   INT   Connectivity of block in donor mesh (1:nelndb,1:numebb)
C  NDLSTB  INT   Vector of nodes in element block donor mesh (1:nodesb)

C     ******************************************************************

      character*(32) typ,typa

      include 'bmesh.blk'
      include 'ebbyeb.blk'
      include 'ex2tp.blk'

      DIMENSION ICONB(NELNDB,*),NDLSTB(*)

C     ******************************************************************

C     READ ELEMENT NAMES AND ID BLOCKS

      CALL EXGELB(NTP3EX,IDBLKB,TYP,NUMEBB,NELNDB,NATRIB,IERR)
      CALL EXUPCS(TYP)
      CALL EXGELB(NTP2EX,IDBLKA,TYPA,IDUM,IDUM,IDUM,IERR)
      CALL EXUPCS(TYPA)

c check here for match of mesh-A element block to mesh-B element block.
c Probably need to put this into a DO LOOP over all element block id's
c in mesh B. For now, assume that element blocks match if id's match.
c Only check element types.

      IF (TYP(1:3) .NE. TYPA(1:3)) THEN
        CALL ERROR('RDB2',
     &    'ELEMENT TYPE MISMATCH - MESH-B DOES NOT MATCH MESH-A',
     &    'ELEMENT BLOCK', IDBLKA,' ',0,
     &    'Execution will continue, but verify that results are OK',
     &    ' ',0)
      END IF

      CALL EXGELC(NTP3EX,IDBLKB,ICONB(1,1),IERR)

      DO 5 I = 1,NODESB
        NDLSTB(I) = 0
 5    CONTINUE

      DO 10 IEL = 1, NUMEBB
        DO 20 INODE = 1, NELNDB
          NDLSTB(ICONB(INODE,IEL)) = 1
 20     CONTINUE
 10   CONTINUE

      NUMNDB = 0

      DO 30 I = 1, NODESB
        IF (NDLSTB(I) .EQ. 1)THEN
          NUMNDB = NUMNDB + 1
          NDLSTB(NUMNDB) = I
        END IF
 30   CONTINUE

      RETURN
      END
