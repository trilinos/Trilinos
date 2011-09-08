C Copyright (c) 1998 Sandia Corporation. Under the terms of Contract
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
C=================================================================
C
C This is just a simple test program to test the fortran interface
C for the NEMESIS I library.
C
C This file was created by translating ne_test.c into fortran.
C
C=================================================================

C=================================================================
      PROGRAM NETEST
C=================================================================
      INCLUDE 'exodusII.inc'
      INCLUDE 'ne_test.inc'

C local variables
      INTEGER NEID, IO_WS, CPU_WS, T_PASS, T_FAIL, DBG_FLAG, IERR
      CHARACTER FNAME*256, YO*6
      REAL VERSION

      YO = 'NETEST'
      IO_WS = 0
      CPU_WS = 0
      T_PASS = 0
      T_FAIL = 0
      DBG_FLAG = 0

C now let's get going...

C I don't care about input arguements, so the file name will be ne_test.nemI
      FNAME = 'ne_test.nemI'
C and set the debug flag to 0
      DBG_FLAG = 0

      PRINT*, '******************Output Tests*****************'
C create the exodus II file
      PRINT*, 'creating ExodusII file...'
      NEID = EXCRE(FNAME, EXCLOB, CPU_WS, IO_WS, IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        PRINT*, YO, ': ERROR, unable to create test file', FNAME, '!'
        GOTO 100
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF

C Test the output of initial information
      PRINT*, 'testing init info output...'
      CALL NEPII(NEID, NPROC, NPROCF, 'S', IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        IF (DBG_FLAG.EQ.1) THEN
          GOTO 100
        END IF
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF

C Test the output of initial global information
      PRINT*, 'testing global init info output...'
      CALL NEPIG(NEID, NNG, NEG, NEBG, NNSG, NSSG, IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        IF (DBG_FLAG.EQ.1) THEN
          GOTO 100
        END IF
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF

C Test the output of the global element block IDs
      PRINT*, 'testing global element block ID output...'
      CALL NETPEBI(NEID, IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        IF (DBG_FLAG.EQ.1) THEN
          GOTO 100
        END IF
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF

C Test the output of the global node-set info
      PRINT*, 'testing global node-set params output...'
      CALL NETPNSP(NEID, IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        IF (DBG_FLAG.EQ.1) THEN
          GOTO 100
        END IF
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF

C Test the output of the global side-set info
      PRINT*, 'testing global side-set params output...'
      CALL NETPSSP(NEID, IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        IF (DBG_FLAG.EQ.1) THEN
          GOTO 100
        END IF
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF

C Test the output of the concatenated load-balance parameters
      PRINT*, 'testing concatenated load balance info output...'
      CALL NETPLBPC(NEID, IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        IF (DBG_FLAG.EQ.1) THEN
          GOTO 100
        END IF
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF

C Test the output of the node map
      PRINT*, 'testing node map output...'
      CALL NETPNM(NEID, IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        IF (DBG_FLAG.EQ.1) THEN
          GOTO 100
        END IF
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF

C Test the output of the element map
      PRINT*, 'testing element map output...'
      CALL NETPEM(NEID, IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        IF (DBG_FLAG.EQ.1) THEN
          GOTO 100
        END IF
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF

C Test the output of the concatenated communication map params
#ifndef Build64  
      PRINT*, 'testing concatenated communication map params output...'
      CALL NETPCMPC(NEID, IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        IF (DBG_FLAG.EQ.1) THEN
          GOTO 100
        END IF
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF
#endif
      
C Test nodal communication map output
      PRINT*, 'testing nodal communication map output...'
      CALL NETPNCM(NEID, IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        IF (DBG_FLAG.EQ.1) THEN
          GOTO 100
        END IF
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF

C Test elemental communication map output
      PRINT*, 'testing elemental communication map output...'
      CALL NETPECM(NEID, IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        IF (DBG_FLAG.EQ.1) THEN
          GOTO 100
        END IF
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF

C Close the ExodusII/Nemesis test file
      PRINT*, 'closing ExodusII file...'
      CALL EXCLOS(NEID, IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        PRINT*, YO, ': ERROR, unable to close test file', FNAME, '!'
          GOTO 100
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF

C=================================================================
C                       INPUT TEST SECTION
C=================================================================

      PRINT*, '******************Input Tests******************'

C Re-open the ExodusII/NemesisI file 
      PRINT*, 'reopening ExodusII file...'
      NEID =  EXOPEN(FNAME, EXREAD, CPU_WS, IO_WS, VERSION, IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        PRINT*, YO, ': ERROR, unable to open test file', FNAME, '!'
          GOTO 100
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF

C Test read of of the initial information
      PRINT*, 'testing init info input...'
      CALL NETGII(NEID, IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        IF (DBG_FLAG.EQ.1) THEN
          GOTO 100
        END IF
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF

C Test read of initial global information
      PRINT*, 'testing global init info input...'
      CALL NETGIG(NEID, IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        IF (DBG_FLAG.EQ.1) THEN
          GOTO 100
        END IF
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF

C Test read of global element block IDs
      PRINT*, 'testing global element block IDs input...'
      CALL NETGEBI(NEID, IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        IF (DBG_FLAG.EQ.1) THEN
          GOTO 100
        END IF
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF

C Test read of global node-set params
      PRINT*, 'testing global node-set params input...'
      CALL NETGNSP(NEID, IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        IF (DBG_FLAG.EQ.1) THEN
          GOTO 100
        END IF
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF

C Test read of global side-set params
      PRINT*, 'testing global side-set params input...'
      CALL NETGSSP(NEID, IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        IF (DBG_FLAG.EQ.1) THEN
          GOTO 100
        END IF
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF

C Test read of load-balance params
      PRINT*, 'testing load-balance params input...'
      CALL NETGLBP(NEID, IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        IF (DBG_FLAG.EQ.1) THEN
          GOTO 100
        END IF
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF

C Test read of the node map
      PRINT*, 'testing node map input...'
      CALL NETGNM(NEID, IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        IF (DBG_FLAG.EQ.1) THEN
          GOTO 100
        END IF
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF

C Test read of the element map
      PRINT*, 'testing element map input...'
      CALL NETGEM(NEID, IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        IF (DBG_FLAG.EQ.1) THEN
          GOTO 100
        END IF
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF

C Test read of nodal communication maps
      PRINT*, 'testing nodal communication map input...'
      CALL NETGNCM(NEID, IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        IF (DBG_FLAG.EQ.1) THEN
          GOTO 100
        END IF
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF

C Test read of elemental communication maps
      PRINT*, 'testing elemental communication map input...'
      CALL NETGECM(NEID, IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        IF (DBG_FLAG.EQ.1) THEN
          GOTO 100
        END IF
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF

C Close the ExodusII/Nemesis test file
      PRINT*, 'closing ExodusII file...'
      CALL EXCLOS(NEID, IERR)
      IF (IERR.NE.0) THEN
        PRINT*, '...FAILED'
        T_FAIL = T_FAIL + 1
        PRINT*, YO, ': ERROR, unable to close test file', FNAME, '!'
          GOTO 100
      ELSE
        PRINT*, '...successful'
        T_PASS = T_PASS + 1
      END IF

      PRINT*, 'Tests Passed: ', T_PASS
      PRINT*, 'Tests Failed: ', T_FAIL

  100 CONTINUE
      END

C=================================================================
      SUBROUTINE NETPEBI(NEID, IERR)
C=================================================================

      INCLUDE 'ne_test.inc'

      INTEGER I, EBLK_IDS(NEBG)
      INTEGER EBLK_CNTS(NEBG)

      DO 110 I=1,NEBG
        EBLK_IDS(I) = I
        EBLK_CNTS(I) = 10
  110 CONTINUE

      CALL NEPEBIG(NEID, EBLK_IDS, EBLK_CNTS, IERR)

      END

C=================================================================
      SUBROUTINE NETPNSP(NEID, IERR)
C=================================================================

      INCLUDE 'ne_test.inc'

      INTEGER I, GLBL_IDS(NNSG), GLBL_NC(NNSG), GLBL_DFC(NNSG)

      DO 120 I = 1,NNSG
        GLBL_IDS(I) = 2 * I
        GLBL_NC(I) = 3 * I
        GLBL_DFC(I) = 1
  120 CONTINUE

      CALL NEPNSPG(NEID, GLBL_IDS, GLBL_NC, GLBL_DFC, IERR)

      END

C=================================================================
      SUBROUTINE NETPSSP(NEID, IERR)
C=================================================================

      INCLUDE 'ne_test.inc'

      INTEGER I, GLBL_IDS(NSSG), GLBL_ELC(NSSG), GLBL_DFC(NSSG)

      DO 130 I = 1,NSSG
        GLBL_IDS(I) = 3 * I
        GLBL_ELC(I) = 2 * I
        GLBL_DFC(I) = 1
  130 CONTINUE

      CALL NEPSSPG(NEID, GLBL_IDS, GLBL_ELC, GLBL_DFC, IERR)

      END

C=================================================================
      SUBROUTINE NETPLBPC(NEID, IERR)
C=================================================================

      INCLUDE 'ne_test.inc'

      INTEGER IPROC, NUM_IN(NPROCF), NUM_BN(NPROCF), NUM_EN(NPROCF),
     1 NUM_IE(NPROCF), NUM_BE(NPROCF), NUM_NCM(NPROCF), NUM_ECM(NPROCF)

      DO 140 IPROC = 1,NPROCF
        NUM_IN(IPROC) = NINTN
        NUM_BN(IPROC) = NBORN
        NUM_EN(IPROC) = NEXTN

        NUM_IE(IPROC) = NINTE
        NUM_BE(IPROC) = NBORE

        NUM_NCM(IPROC) = NNCMAP
        NUM_ECM(IPROC) = NECMAP
  140 CONTINUE

      CALL NEPLBPC(NEID, NUM_IN, NUM_BN, NUM_EN, NUM_IE, NUM_BE,
     1 NUM_NCM, NUM_ECM, IERR)

      END

C=================================================================
      SUBROUTINE NETPNM(NEID, IERR)
C=================================================================

      INCLUDE 'ne_test.inc'

      INTEGER IPROC, I, J, NMAPI(NINTN), NMAPB(NBORN), NMAPE(NEXTN)

      I = 0
      DO 200 IPROC = 0,(NPROCF-1)
        DO 150 J = 1,NINTN
          NMAPI(J) = I
          I = I + 1
  150   CONTINUE
        DO 160 J = 1,NBORN
          NMAPB(J) = I
          I = I + 1
  160   CONTINUE
        DO 170 J = 1,NEXTN
          NMAPE(J) = I
          I = I + 1
  170   CONTINUE

        I = 0

        CALL NEPNM(NEID, NMAPI, NMAPB, NMAPE, IPROC, IERR)
        IF (IERR.NE.0) GOTO 210

  200 CONTINUE

  210 CONTINUE
      END

C=================================================================
      SUBROUTINE NETPEM(NEID, IERR)
C=================================================================

      INCLUDE 'ne_test.inc'

      INTEGER IPROC, I, J, EMAPI(NINTE), EMAPB(NBORE)

      I = 0
      DO 200 IPROC = 0,(NPROCF-1)
        DO 150 J = 1,NINTE
          EMAPI(J) = I
          I = I + 1
  150   CONTINUE
        DO 160 J = 1,NBORE
          EMAPB(J) = I
          I = I + 1
  160   CONTINUE

        I = 0

        CALL NEPEM(NEID, EMAPI, EMAPB, IPROC, IERR)
        IF (IERR.NE.0) GOTO 210

  200 CONTINUE

  210 CONTINUE
      END

C=================================================================
      SUBROUTINE NETPCMPC(NEID, IERR)
C=================================================================

      INCLUDE 'ne_test.inc'

      INTEGER IPROC, I, NCNTR, ECNTR, NMAPIDS(NNCXNPF),
     1 NMAPCNT(NNCXNPF), NMAPPROC(NPROCF+1), EMAPIDS(NECXNPF),
     1 EMAPCNT(NECXNPF), EMAPPROC(NPROCF+1)

      NMAPPROC(1) = 0
      EMAPPROC(1) = 0
      NCNTR = 1
      ECNTR = 1
      DO 200 IPROC = 1,NPROCF
        DO 150 I = 1,NNCMAP
          NMAPIDS(NCNTR) = I
          NMAPCNT(NCNTR) = NCNTCM
          NCNTR = NCNTR + 1
  150   CONTINUE
        DO 160 I = 1,NECMAP
          EMAPIDS(ECNTR) = 2*I
          EMAPCNT(ECNTR) = ECNTCM
          ECNTR = ECNTR + 1
  160   CONTINUE

      NMAPPROC(IPROC+1) = NMAPPROC(IPROC) + NNCMAP
      EMAPPROC(IPROC+1) = EMAPPROC(IPROC) + NECMAP

  200 CONTINUE

      CALL NEPCMPC(NEID, NMAPIDS, NMAPCNT, NMAPPROC, EMAPIDS, EMAPCNT, 
     1 EMAPPROC, IERR)

      END

C=================================================================
      SUBROUTINE NETPNCM(NEID, IERR)
C=================================================================

      INCLUDE 'ne_test.inc'

      INTEGER IPROC, I, NMAPIDS(NNCMAP), NIDS(NCNTCM), PIDS(NCNTCM)

      DO 200 IPROC = 0,(NPROCF-1)
        DO 150 I = 1,NNCMAP
          NMAPIDS(I) = I
  150   CONTINUE
        DO 160 I = 1,NCNTCM
          NIDS(I) = 2*I
          PIDS(I) = 3*I
  160   CONTINUE

        DO 170 I=1,NNCMAP
          CALL NEPNCM(NEID, NMAPIDS(I), NIDS, PIDS, IPROC, IERR)
          IF (IERR.NE.0) GOTO 210
  170   CONTINUE

  200 CONTINUE

  210 CONTINUE
      END

C=================================================================
      SUBROUTINE NETPECM(NEID, IERR)
C=================================================================

      INCLUDE 'ne_test.inc'

      INTEGER IPROC, I, EMAPIDS(NECMAP), EIDS(ECNTCM), PIDS(ECNTCM),
     1 SIDS(ECNTCM)

      DO 200 IPROC = 0,(NPROCF-1)
        DO 150 I = 1,NECMAP
          EMAPIDS(I) = 2*I
  150   CONTINUE
        DO 160 I = 1,ECNTCM
          EIDS(I) = 2*I
          SIDS(I) = 3*I
          PIDS(I) = 4*I
  160   CONTINUE

        DO 170 I=1,NECMAP
          CALL NEPECM(NEID, EMAPIDS(I), EIDS, SIDS, PIDS, IPROC, IERR)
          IF (IERR.NE.0) GOTO 210
  170   CONTINUE

  200 CONTINUE

  210 CONTINUE
      END

C=================================================================
      SUBROUTINE NETGII(NEID, IERR)
C=================================================================

      INCLUDE 'ne_test.inc'

      INTEGER NP, NPF
      CHARACTER FTYPE*2

      CALL NEGII(NEID, NP, NPF, FTYPE, IERR)

      IF (IERR.NE.0) GOTO 210

      IF (NP.NE.NPROC) IERR = -1
      IF (NPF.NE.NPROCF) IERR = -1
      IF (NP.NE.NPROC) IERR = -1

  210 CONTINUE
      END

C=================================================================
      SUBROUTINE NETGIG(NEID, IERR)
C=================================================================

      INCLUDE 'ne_test.inc'

      INTEGER NUMNG, NUMEG, NUMEBG, NUMNSG, NUMSSG

      CALL NEGIG(NEID, NUMNG, NUMEG, NUMEBG, NUMNSG, NUMSSG, IERR)

      IF (IERR.NE.0) GOTO 210

      IF (NUMNG.NE.NNG) IERR = -1
      IF (NUMEG.NE.NEG) IERR = -1
      IF (NUMEBG.NE.NEBG) IERR = -1
      IF (NUMNSG.NE.NNSG) IERR = -1
      IF (NUMSSG.NE.NSSG) IERR = -1

  210 CONTINUE
      END

C=================================================================
      SUBROUTINE NETGEBI(NEID, IERR)
C=================================================================

      INCLUDE 'ne_test.inc'

      INTEGER I, EBLK_IDS(NEBG)
      INTEGER EBLK_CNTS(NEBG)

      CALL NEGEBIG(NEID, EBLK_IDS, EBLK_CNTS, IERR)

      IF (IERR.NE.0) GOTO 210

      DO 150 I=1,NEBG
        IF (EBLK_IDS(I).NE.I) IERR = -1
        IF (EBLK_CNTS(I) .NE. 10) IERR = -1
  150 CONTINUE

  210 CONTINUE
      END

C=================================================================
      SUBROUTINE NETGNSP(NEID, IERR)
C=================================================================

      INCLUDE 'ne_test.inc'

      INTEGER I, GLBL_IDS(NNSG), GLBL_NC(NNSG), GLBL_DFC(NNSG)

      CALL NEGNSPG(NEID, GLBL_IDS, GLBL_NC, GLBL_DFC, IERR)

      IF (IERR.NE.0) GOTO 210

      DO 150 I=1,NNSG
        IF (GLBL_IDS(I).NE.(2*I)) IERR = -1
        IF (GLBL_NC(I).NE.(3*I)) IERR = -1
        IF (GLBL_DFC(I).NE.1) IERR = -1
  150 CONTINUE

  210 CONTINUE
      END

C=================================================================
      SUBROUTINE NETGSSP(NEID, IERR)
C=================================================================

      INCLUDE 'ne_test.inc'

      INTEGER I, GLBL_IDS(NSSG), GLBL_EC(NSSG), GLBL_DFC(NSSG)

      CALL NEGSSPG(NEID, GLBL_IDS, GLBL_EC, GLBL_DFC, IERR)

      IF (IERR.NE.0) GOTO 210

      DO 150 I=1,NNSG
        IF (GLBL_IDS(I).NE.(3*I)) IERR = -1
        IF (GLBL_EC(I).NE.(2*I)) IERR = -1
        IF (GLBL_DFC(I).NE.1) IERR = -1
  150 CONTINUE

  210 CONTINUE
      END

C=================================================================
      SUBROUTINE NETGLBP(NEID, IERR)
C=================================================================

      INCLUDE 'ne_test.inc'

      INTEGER IPROC, NUM_IN, NUM_BN, NUM_EN, NUM_IE, NUM_BE,
     *  NUM_NCM, NUM_ECM

      DO 150 IPROC = 0,(NPROCF-1)
        CALL NEGLBP(NEID, NUM_IN, NUM_BN, NUM_EN, NUM_IE, NUM_BE,
     1   NUM_NCM, NUM_ECM, IPROC, IERR)

      IF (IERR.NE.0) GOTO 210

        IF(NUM_IN.NE.NINTN) IERR = -1
        IF(NUM_BN.NE.NBORN) IERR = -1
        IF(NUM_EN.NE.NEXTN) IERR = -1
        IF(NUM_IE.NE.NINTE) IERR = -1
        IF(NUM_BE.NE.NBORE) IERR = -1
        IF(NUM_NCM.NE.NNCMAP) IERR = -1
        IF(NUM_ECM.NE.NECMAP) IERR = -1
  150 CONTINUE

  210 CONTINUE
      END

C=================================================================
      SUBROUTINE NETGNM(NEID, IERR)
C=================================================================

      INCLUDE 'ne_test.inc'

      INTEGER IPROC, I, J, NMAPI(NINTN), NMAPB(NBORN), NMAPE(NEXTN)

      I = 0
      DO 200 IPROC = 0,(NPROCF-1)

        CALL NEGNM(NEID, NMAPI, NMAPB, NMAPE, IPROC, IERR)

        IF (IERR.NE.0) GOTO 210

        DO 150 J = 1,NINTN
          IF (NMAPI(J).NE.I) ERR = -1
          I = I + 1
  150   CONTINUE
        DO 160 J = 1,NBORN
          IF (NMAPB(J).NE.I) ERR = -1
          I = I + 1
  160   CONTINUE
        DO 170 J = 1,NEXTN
          IF (NMAPE(J).NE.I) ERR = -1
          I = I + 1
  170   CONTINUE

        I = 0

        IF (IERR.NE.0) GOTO 210

  200 CONTINUE

  210 CONTINUE
      END

C=================================================================
      SUBROUTINE NETGEM(NEID, IERR)
C=================================================================

      INCLUDE 'ne_test.inc'

      INTEGER IPROC, I, J, EMAPI(NINTE), EMAPB(NBORE)

      I = 0
      DO 200 IPROC = 0,(NPROCF-1)
        CALL NEGEM(NEID, EMAPI, EMAPB, IPROC, IERR)

        IF (IERR.NE.0) GOTO 210

        DO 150 J = 1,NINTE
          IF (EMAPI(J).NE.I) ERR = -1
          I = I + 1
  150   CONTINUE
        DO 160 J = 1,NBORE
          IF (EMAPB(J).NE.I) ERR = -1
          I = I + 1
  160   CONTINUE

        I = 0

        IF (IERR.NE.0) GOTO 210

  200 CONTINUE

  210 CONTINUE
      END

C=================================================================
      SUBROUTINE NETGNCM(NEID, IERR)
C=================================================================

      INCLUDE 'ne_test.inc'

      INTEGER IPROC, I, J, NMAPIDS(NNCMAP), NMAPCNT(NNCMAP),
     1 NIDS(NCNTCM), PIDS(NCNTCM), EMAPIDS(NECMAP), EMAPCNT(NECMAP)

      DO 200 IPROC = 0,(NPROCF-1)
        CALL NEGCMP(NEID, NMAPIDS, NMAPCNT, EMAPIDS, EMAPCNT,
     1   IPROC, IERR)

        IF (IERR.NE.0) GOTO 210

        DO 170 I = 1,NNCMAP
          CALL NEGNCM(NEID, NMAPIDS(I), NIDS, PIDS, IPROC, IERR)

          IF (IERR.NE.0) GOTO 210

          IF (NMAPIDS(I).NE.I) IERR = -1
          DO 160 J = 1,NCNTCM
            IF (NIDS(J).NE.2*J) IERR = -1
            IF (PIDS(J).NE.3*J) IERR = -1
  160     CONTINUE

          IF (IERR.NE.0) GOTO 210
  170   CONTINUE

  200 CONTINUE

  210 CONTINUE
      END

C=================================================================
      SUBROUTINE NETGECM(NEID, IERR)
C=================================================================

      INCLUDE 'ne_test.inc'

      INTEGER IPROC, I, EMAPIDS(NECMAP), EMAPCNT(NECMAP), EIDS(ECNTCM),
     1 PIDS(ECNTCM), SIDS(ECNTCM), NMAPIDS(NNCMAP), NMAPCNT(NNCMAP)

      DO 200 IPROC = 0,(NPROCF-1)
        CALL NEGCMP(NEID, NMAPIDS, NMAPCNT, EMAPIDS, EMAPCNT,
     1   IPROC, IERR)

        IF (IERR.NE.0) GOTO 210

        DO 170 I = 1,NECMAP
          CALL NEGECM(NEID, EMAPIDS(I), EIDS, SIDS, PIDS, IPROC, IERR)

          IF (IERR.NE.0) GOTO 210

          IF (EMAPIDS(I).NE.(2*I)) IERR = -1
          DO 160 J = 1,ECNTCM
            IF (EIDS(J).NE.2*J) IERR = -1
            IF (SIDS(J).NE.3*J) IERR = -1
            IF (PIDS(J).NE.4*J) IERR = -1
  160     CONTINUE

          IF (IERR.NE.0) GOTO 210
  170   CONTINUE

  200 CONTINUE

  210 CONTINUE
      END

