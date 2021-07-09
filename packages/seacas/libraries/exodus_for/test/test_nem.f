C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details
C=================================================================

C This is just a simple test program to test the fortran interface
C for the NEMESIS I library.

C This file was created by translating ne_test.c into fortran.

C=================================================================

C=================================================================
      PROGRAM NETEST
C=================================================================
      INCLUDE 'exodusII.inc'
      INCLUDE 'test_nem.inc'

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

C I don't care about input arguments, so the file name will be ne_test.nemI
      FNAME = 'test_nem.exo'
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
      CALL EXPII(NEID, NPROC, NPROCF, 'S', IERR)
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
      CALL EXPIG(NEID, NNG, NEG, NEBG, NNSG, NSSG, IERR)
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
      CALL EXTPEBI(NEID, IERR)
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
      CALL EXTPNSP(NEID, IERR)
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
      CALL EXTPSSP(NEID, IERR)
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
      CALL EXTPLBPC(NEID, IERR)
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
      CALL EXTPNM(NEID, IERR)
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
      CALL EXTPEM(NEID, IERR)
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
      PRINT*, 'testing concatenated communication map params output...'
      CALL EXTPCMPC(NEID, IERR)
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

C Test nodal communication map output
      PRINT*, 'testing nodal communication map output...'
      CALL EXTPNCM(NEID, IERR)
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
      CALL EXTPECM(NEID, IERR)
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
      CALL EXTGII(NEID, IERR)
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
      CALL EXTGIG(NEID, IERR)
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
      CALL EXTGEBI(NEID, IERR)
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
      CALL EXTGNSP(NEID, IERR)
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
      CALL EXTGSSP(NEID, IERR)
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
      CALL EXTGLBP(NEID, IERR)
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
      CALL EXTGNM(NEID, IERR)
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
      CALL EXTGEM(NEID, IERR)
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
      CALL EXTGNCM(NEID, IERR)
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
      CALL EXTGECM(NEID, IERR)
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
      SUBROUTINE EXTPEBI(NEID, IERR)
C=================================================================

      INCLUDE 'test_nem.inc'

      INTEGER I, EBLK_IDS(NEBG)
      INTEGER EBLK_CNTS(NEBG)

      DO 110 I=1,NEBG
        EBLK_IDS(I) = I
        EBLK_CNTS(I) = 10
  110 CONTINUE

      CALL EXPEBIG(NEID, EBLK_IDS, EBLK_CNTS, IERR)

      END

C=================================================================
      SUBROUTINE EXTPNSP(NEID, IERR)
C=================================================================

      INCLUDE 'test_nem.inc'

      INTEGER I, GLBL_IDS(NNSG), GLBL_NC(NNSG), GLBL_DFC(NNSG)

      DO 120 I = 1,NNSG
        GLBL_IDS(I) = 2 * I
        GLBL_NC(I) = 3 * I
        GLBL_DFC(I) = 1
  120 CONTINUE

      CALL EXPNSPG(NEID, GLBL_IDS, GLBL_NC, GLBL_DFC, IERR)

      END

C=================================================================
      SUBROUTINE EXTPSSP(NEID, IERR)
C=================================================================

      INCLUDE 'test_nem.inc'

      INTEGER I, GLBL_IDS(NSSG), GLBL_ELC(NSSG), GLBL_DFC(NSSG)

      DO 130 I = 1,NSSG
        GLBL_IDS(I) = 3 * I
        GLBL_ELC(I) = 2 * I
        GLBL_DFC(I) = 1
  130 CONTINUE

      CALL EXPSSPG(NEID, GLBL_IDS, GLBL_ELC, GLBL_DFC, IERR)

      END

C=================================================================
      SUBROUTINE EXTPLBPC(NEID, IERR)
C=================================================================

      INCLUDE 'test_nem.inc'

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

      CALL EXPLBPC(NEID, NUM_IN, NUM_BN, NUM_EN, NUM_IE, NUM_BE,
     1 NUM_NCM, NUM_ECM, IERR)

      END

C=================================================================
      SUBROUTINE EXTPNM(NEID, IERR)
C=================================================================

      INCLUDE 'test_nem.inc'

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

        CALL EXPNMP(NEID, NMAPI, NMAPB, NMAPE, IPROC, IERR)
        IF (IERR.NE.0) GOTO 210

  200 CONTINUE

  210 CONTINUE
      END

C=================================================================
      SUBROUTINE EXTPEM(NEID, IERR)
C=================================================================

      INCLUDE 'test_nem.inc'

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

        CALL EXPEMP(NEID, EMAPI, EMAPB, IPROC, IERR)
        IF (IERR.NE.0) GOTO 210

  200 CONTINUE

  210 CONTINUE
      END

C=================================================================
      SUBROUTINE EXTPCMPC(NEID, IERR)
C=================================================================

      INCLUDE 'test_nem.inc'

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

      CALL EXPCMPC(NEID, NMAPIDS, NMAPCNT, NMAPPROC, EMAPIDS, EMAPCNT,
     1 EMAPPROC, IERR)

      END

C=================================================================
      SUBROUTINE EXTPNCM(NEID, IERR)
C=================================================================

      INCLUDE 'test_nem.inc'

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
          CALL EXPNCM(NEID, NMAPIDS(I), NIDS, PIDS, IPROC, IERR)
          IF (IERR.NE.0) GOTO 210
  170   CONTINUE

  200 CONTINUE

  210 CONTINUE
      END

C=================================================================
      SUBROUTINE EXTPECM(NEID, IERR)
C=================================================================

      INCLUDE 'test_nem.inc'

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
          CALL EXPECM(NEID, EMAPIDS(I), EIDS, SIDS, PIDS, IPROC, IERR)
          IF (IERR.NE.0) GOTO 210
  170   CONTINUE

  200 CONTINUE

  210 CONTINUE
      END

C=================================================================
      SUBROUTINE EXTGII(NEID, IERR)
C=================================================================

      INCLUDE 'test_nem.inc'

      INTEGER NP, NPF
      CHARACTER FTYPE*2

      CALL EXGII(NEID, NP, NPF, FTYPE, IERR)

      IF (IERR.NE.0) GOTO 210

      IF (NP.NE.NPROC) IERR = -1
      IF (NPF.NE.NPROCF) IERR = -1
      IF (NP.NE.NPROC) IERR = -1

  210 CONTINUE
      END

C=================================================================
      SUBROUTINE EXTGIG(NEID, IERR)
C=================================================================

      INCLUDE 'test_nem.inc'

      INTEGER NUMNG, NUMEG, NUMEBG, NUMNSG, NUMSSG

      CALL EXGIG(NEID, NUMNG, NUMEG, NUMEBG, NUMNSG, NUMSSG, IERR)

      IF (IERR.NE.0) GOTO 210

      IF (NUMNG.NE.NNG) IERR = -1
      IF (NUMEG.NE.NEG) IERR = -1
      IF (NUMEBG.NE.NEBG) IERR = -1
      IF (NUMNSG.NE.NNSG) IERR = -1
      IF (NUMSSG.NE.NSSG) IERR = -1

  210 CONTINUE
      END

C=================================================================
      SUBROUTINE EXTGEBI(NEID, IERR)
C=================================================================

      INCLUDE 'test_nem.inc'

      INTEGER I, EBLK_IDS(NEBG)
      INTEGER EBLK_CNTS(NEBG)

      CALL EXGEBIG(NEID, EBLK_IDS, EBLK_CNTS, IERR)

      IF (IERR.NE.0) GOTO 210

      DO 150 I=1,NEBG
        IF (EBLK_IDS(I).NE.I) IERR = -1
        IF (EBLK_CNTS(I) .NE. 10) IERR = -1
  150 CONTINUE

  210 CONTINUE
      END

C=================================================================
      SUBROUTINE EXTGNSP(NEID, IERR)
C=================================================================

      INCLUDE 'test_nem.inc'

      INTEGER I, GLBL_IDS(NNSG), GLBL_NC(NNSG), GLBL_DFC(NNSG)

      CALL EXGNSPG(NEID, GLBL_IDS, GLBL_NC, GLBL_DFC, IERR)

      IF (IERR.NE.0) GOTO 210

      DO 150 I=1,NNSG
        IF (GLBL_IDS(I).NE.(2*I)) IERR = -1
        IF (GLBL_NC(I).NE.(3*I)) IERR = -1
        IF (GLBL_DFC(I).NE.1) IERR = -1
  150 CONTINUE

  210 CONTINUE
      END

C=================================================================
      SUBROUTINE EXTGSSP(NEID, IERR)
C=================================================================

      INCLUDE 'test_nem.inc'

      INTEGER I, GLBL_IDS(NSSG), GLBL_EC(NSSG), GLBL_DFC(NSSG)

      CALL EXGSSPG(NEID, GLBL_IDS, GLBL_EC, GLBL_DFC, IERR)

      IF (IERR.NE.0) GOTO 210

      DO 150 I=1,NNSG
        IF (GLBL_IDS(I).NE.(3*I)) IERR = -1
        IF (GLBL_EC(I).NE.(2*I)) IERR = -1
        IF (GLBL_DFC(I).NE.1) IERR = -1
  150 CONTINUE

  210 CONTINUE
      END

C=================================================================
      SUBROUTINE EXTGLBP(NEID, IERR)
C=================================================================

      INCLUDE 'test_nem.inc'

      INTEGER IPROC, NUM_IN, NUM_BN, NUM_EN, NUM_IE, NUM_BE,
     *  NUM_NCM, NUM_ECM

      DO 150 IPROC = 0,(NPROCF-1)
        CALL EXGLBP(NEID, NUM_IN, NUM_BN, NUM_EN, NUM_IE, NUM_BE,
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
      SUBROUTINE EXTGNM(NEID, IERR)
C=================================================================

      INCLUDE 'test_nem.inc'

      INTEGER IPROC, I, J, NMAPI(NINTN), NMAPB(NBORN), NMAPE(NEXTN)

      I = 0
      DO 200 IPROC = 0,(NPROCF-1)

        CALL EXGNMP(NEID, NMAPI, NMAPB, NMAPE, IPROC, IERR)

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
      SUBROUTINE EXTGEM(NEID, IERR)
C=================================================================

      INCLUDE 'test_nem.inc'

      INTEGER IPROC, I, J, EMAPI(NINTE), EMAPB(NBORE)

      I = 0
      DO 200 IPROC = 0,(NPROCF-1)
        CALL EXGEMP(NEID, EMAPI, EMAPB, IPROC, IERR)

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
      SUBROUTINE EXTGNCM(NEID, IERR)
C=================================================================

      INCLUDE 'test_nem.inc'

      INTEGER IPROC, I, J, NMAPIDS(NNCMAP), NMAPCNT(NNCMAP),
     1 NIDS(NCNTCM), PIDS(NCNTCM), EMAPIDS(NECMAP), EMAPCNT(NECMAP)

      DO 200 IPROC = 0,(NPROCF-1)
        CALL EXGCMP(NEID, NMAPIDS, NMAPCNT, EMAPIDS, EMAPCNT,
     1   IPROC, IERR)

        IF (IERR.NE.0) GOTO 210

        DO 170 I = 1,NNCMAP
          CALL EXGNCM(NEID, NMAPIDS(I), NIDS, PIDS, IPROC, IERR)

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
      SUBROUTINE EXTGECM(NEID, IERR)
C=================================================================

      INCLUDE 'test_nem.inc'

      INTEGER IPROC, I, EMAPIDS(NECMAP), EMAPCNT(NECMAP), EIDS(ECNTCM),
     1 PIDS(ECNTCM), SIDS(ECNTCM), NMAPIDS(NNCMAP), NMAPCNT(NNCMAP)

      DO 200 IPROC = 0,(NPROCF-1)
        CALL EXGCMP(NEID, NMAPIDS, NMAPCNT, EMAPIDS, EMAPCNT,
     1   IPROC, IERR)

        IF (IERR.NE.0) GOTO 210

        DO 170 I = 1,NECMAP
          CALL EXGECM(NEID, EMAPIDS(I), EIDS, SIDS, PIDS, IPROC, IERR)

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

