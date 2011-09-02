C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.  
C 
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C 

C=======================================================================
*DECK,RDB2
      SUBROUTINE RDB2 (IDBLKB,IDBLKA,ICONB,NDLSTB)
C
C     ******************************************************************
C
C     SUBROUTINE TO GET MESH B DATA
c
C     Calls subroutine ERROR
C
C     Called by MAPVAR
C
C     ******************************************************************
C
C  IDBLKA  INT   Element block I.D. - donor mesh
C  IDBLKB  INT   Element block I.D. - recipient mesh
C  ICONB   INT   Connectivity of block in donor mesh (1:nelndb,1:numebb)
C  NDLSTB  INT   Vector of nodes in element block donor mesh (1:nodesb)
C
C     ******************************************************************
C
      include 'exodusII.inc'
      character*(MXSTLN) typ,typa
C
      include 'bmesh.blk'
      include 'ebbyeb.blk'
      include 'ex2tp.blk'
C
      DIMENSION ICONB(NELNDB,*),NDLSTB(*)
C
C     ******************************************************************
C
C     READ ELEMENT NAMES AND ID BLOCKS
C
      CALL EXGELB(NTP3EX,IDBLKB,TYP,NUMEBB,NELNDB,NATRIB,IERR)
      CALL EXUPCS(TYP)
      CALL EXGELB(NTP2EX,IDBLKA,TYPA,IDUM,IDUM,IDUM,IERR)
      CALL EXUPCS(TYPA)
c
c check here for match of mesh-A element block to mesh-B element block. 
c Probably need to put this into a DO LOOP over all element block id's
c in mesh B. For now, assume that element blocks match if id's match. 
c Only check element types.
c
      IF (TYP(1:3) .NE. TYPA(1:3))CALL ERROR('RDB2','ELEMENT TYPE 
     &  MISMATCH - MESH-B DOES NOT MATCH MESH-A','ELEMENT BLOCK',
     &  IDBLKA,' ',0,' ',' ',1)
C
C
      CALL EXGELC(NTP3EX,IDBLKB,ICONB(1,1),IERR)
C
      DO 5 I = 1,NODESB
        NDLSTB(I) = 0
 5    CONTINUE
C
      DO 10 IEL = 1, NUMEBB
        DO 20 INODE = 1, NELNDB
          NDLSTB(ICONB(INODE,IEL)) = 1
 20     CONTINUE
 10   CONTINUE
C
      NUMNDB = 0
C
      DO 30 I = 1, NODESB
        IF (NDLSTB(I) .EQ. 1)THEN
          NUMNDB = NUMNDB + 1
          NDLSTB(NUMNDB) = I
        END IF
 30   CONTINUE
c
      RETURN
      END
