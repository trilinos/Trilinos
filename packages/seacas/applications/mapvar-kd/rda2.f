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
*DECK,RDA2
      SUBROUTINE RDA2 (IDBLKA,ICONA,NDLSTA,STATUS,MAXLN)
C
C     ******************************************************************
C
C     SUBROUTINE TO READ MESH A, WRITE MESH C DATA AS APPROPRIATE
C
C     Calls subroutine ERROR
C
C     Called by MAPVAR
C
C     ******************************************************************
C
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
C
C     ******************************************************************
C
      include 'exodusII.inc'
      CHARACTER*(MXSTLN) TYP
C
      include 'aexds1.blk'
      include 'aexds2.blk'
      include 'amesh.blk'
      include 'ebbyeb.blk'
      include 'ex2tp.blk'
      include 'steps.blk'
      include 'tapes.blk'
      include 'debg.blk'
C
      DIMENSION ICONA(NELNDA,*),NDLSTA(*),STATUS(*)
C
C     ******************************************************************
C
C element type per element block
C
C
C fix this routine when i have time
C create array 1-nnodes
C loop over all elements - add 1 to value in array whenever 
C                          node appears in connectivity 
C maxln = max value of array
C
        CALL EXGELB (NTP2EX,IDBLKA,TYP,NUMEBA,NELNDA,
     &               NATRIB,IERR)
        CALL EXUPCS(TYP)
c
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
C
C
      CALL EXGELC(NTP2EX,IDBLKA,ICONA(1,1),IERR)
C
      DO 5 I = 1, NODESA
        NDLSTA(I) = 0
 5    CONTINUE
C
      DO 10 IEL = 1, NUMEBA
        DO 20 INODE = 1, NELNDA
          NDLSTA(ICONA(INODE,IEL)) = NDLSTA(ICONA(INODE,IEL)) + 1
 20     CONTINUE
 10   CONTINUE
C
      NUMNDA = 0
C
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
      if (idebug .ge. 1) then
        write(nout, *)
     $        '        In RDA2 -- Max elements per node = ', maxln
        write(ntpout, *)
     $       '        In RDA2 -- Max elements per node = ', maxln
      end if
      
C
C get STATUS array for use in SEARCH so that dead elements can be 
C eliminated from the search
C
      DO 99 I = 1, NUMEBA
         STATUS(I) = 0.
   99 CONTINUE
      DO 100 ISTATUS = 1, NVAREL
         IF (NAMVAR(nvargp+ISTATUS) .NE. 'STATUS')GO TO 100
         CALL EXGEV(NTP2EX,ISTEP,ISTATUS,IDBLKA,NUMEBA,STATUS,IERR)
  100 CONTINUE
c
      RETURN
      END
