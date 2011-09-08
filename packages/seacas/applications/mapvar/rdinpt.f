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
*DECK,RDINPT
      SUBROUTINE RDINPT (TIMES,IDA,IDB,MP,IMP)
C
C     ******************************************************************
C
C     SUBROUTINE TO READ, CHECK AND PRINT INPUT DATA FROM STD-INPUT
C     BATCH TYPE EXECUTION IS ACCOMPLISHEDBY PIPING INPUT DATA
C     FROM A TEXT FILE; NORMALLY INPUT READ INTERACTIVELY.
C     INPUT IS READ UNDER A FREE FIELD FORMAT IN SUBROUTINE FREFLD
C
C     SUBROUTINE FREFLD IS PART OF THE EXTERNAL "SUPES"
C     LIBRARY (SAND86-0911)
C
C     Calls subroutines BANNER, ERROR
C
C     Called by MAPVAR
C
C     ******************************************************************
C
      include 'exodusII.inc'
      CHARACTER*10 CVALUE
C
      include 'amesh.blk'
      include 'bmesh.blk'
      include 'contrl.blk'
      include 'ex2tp.blk'
      include 'rundat.blk'
      include 'schdat.blk'
      include 'steps.blk'
      include 'tapes.blk'
C
      DIMENSION KVALUE(6),CVALUE(6),IVALUE(6),RVALUE(6)
      DIMENSION TIMES(*),IDA(*),IDB(*),MP(3,*)
C
C     ******************************************************************
      MFIELD = 6
C
C     PRINT RUN-TIME DATA
C
      WRITE (NOUT, 1000)
      WRITE (NTPOUT, 1000)
      CALL BANNR2 (84,QAINFO(1),NOUT)
      CALL BANNR2 (84,QAINFO(1),NTPOUT)
      WRITE (NOUT, 1010)
      WRITE (NTPOUT, 1010)
      WRITE (NOUT, 1020) QAINFO(3)
      WRITE (NTPOUT, 1020) QAINFO(3)
      WRITE (NOUT, 1030) QAINFO(5)
      WRITE (NTPOUT, 1030) QAINFO(5)
      WRITE (NOUT, 1040) QAINFO(6)
      WRITE (NTPOUT, 1040) QAINFO(6)
      WRITE (NOUT, 1050)
      WRITE (NTPOUT, 1050)
C
C default map
C
      CALL EXGEBI(NTP2EX,IDA,IERR)
      CALL EXGEBI(NTP3EX,IDB,IERR)
      IMP = 0
      IMP2 = 0
      DO 2 I = 1, NBLKSB
        DO 3 J = 1, NBLKSA
          IF (IDB(I) .EQ. IDA(J))THEN
            IMP = IMP + 1
            MP(1,IMP) = IDA(J)
            MP(2,IMP) = IDB(I)
            MP(3,IMP) = ISCHEM
            GO TO 2
          END IF
    3   CONTINUE
    2 CONTINUE
C
      CALL EXINQ (NTP2EX,EXTIMS,NUMTIM,FDUM,CDUM,IERR)
      CALL EXGATM (NTP2EX,TIMES,IERR)
      OUTTIM = -1.
      ISTEP = NUMTIM
      RTIME = TIMES(NUMTIM)
    4 WRITE (NOUT, 1060)
      WRITE (NOUT, 1061)
C
    5 CALL FREFLD (0,0,'CMD >',MFIELD,IOSTAT,NFIELD,KVALUE,CVALUE,     
     1IVALUE,RVALUE)
C
      IF (KVALUE(1)      .NE. 0)     GO TO 10
      IF (CVALUE(1)(1:3) .EQ. 'HEL') GO TO 20
      IF (CVALUE(1)(1:3) .EQ. 'TIM') GO TO 30
      IF (CVALUE(1)(1:3) .EQ. 'STE') GO TO 35
      IF (CVALUE(1)(1:3) .EQ. 'OUT') GO TO 39
      IF (CVALUE(1)(1:3) .EQ. 'LIS') GO TO 40
      IF (CVALUE(1)(1:3) .EQ. 'SCH') GO TO 50
      IF (CVALUE(1)(1:3) .EQ. 'SEA') GO TO 60
      IF (CVALUE(1)(1:3) .EQ. 'DEF') GO TO 70
      IF (CVALUE(1)(1:3) .EQ. 'MAP') GO TO 80
      IF (CVALUE(1)(1:3) .EQ. 'CHE') GO TO 89
      IF (CVALUE(1)(1:3) .EQ. 'STO') GO TO 90
      IF (CVALUE(1)(1:3) .EQ. 'END') GO TO 100
      IF (CVALUE(1)(1:3) .EQ. 'QUI') GO TO 100
      IF (CVALUE(1)(1:3) .EQ. 'EXI') GO TO 100
      IF (CVALUE(1)(1:3) .EQ. 'RUN') GO TO 100
C
   10 CONTINUE
C
C Bad input
C
      WRITE (NOUT,1100) CVALUE(1)
      WRITE (NTPOUT,1100) CVALUE(1)
      GO TO 5
C
   20 CONTINUE
C
C Help
C
      IF (NFIELD .EQ. 1)THEN
        GO TO 4
      ELSE IF (CVALUE(2)(1:3) .EQ. 'TIM') THEN 
        WRITE(NOUT,2000)
        GO TO 5
      ELSE IF (CVALUE(2)(1:3) .EQ. 'LIS') THEN
        WRITE(NOUT,2010)
        GO TO 5
      ELSE IF (CVALUE(2)(1:3) .EQ. 'SCH') THEN
        WRITE(NOUT,2020)
        GO TO 5
      ELSE IF (CVALUE(2)(1:3) .EQ. 'DEF')THEN
        WRITE(NOUT,2030)
        WRITE(NOUT,2040)
        WRITE(NOUT,2045)
        GO TO 5
      ELSE IF (CVALUE(2)(1:3) .EQ. 'SEA')THEN
        WRITE(NOUT,2050)
        GO TO 5
      ELSE IF (CVALUE(2)(1:3) .EQ. 'MAP')THEN
        WRITE(NOUT,2060)
        GO TO 5
      ELSE IF (CVALUE(2)(1:3) .EQ. 'CHE')THEN
        WRITE(NOUT,2070)
        GO TO 5
      END IF
C
   30 CONTINUE
C
C Time
C
      IF (KVALUE(2) .NE. 1 .AND. KVALUE(2) .NE. 2)THEN
        IF(CVALUE(2)(1:3) .EQ. 'ALL')THEN
          ISTEP = -1
          IDEF = 0
          WRITE(NOUT,3000)
          GO TO 40
        ELSE
          WRITE(NOUT,3010)CVALUE(2)
          GO TO 5
        END IF
      ELSE
C
C convert time to closest time step
C
        RTIME = RVALUE(2)
        ISTEP = NUMTIM
        DO 32 I = 1, NUMTIM - 1
          TMID = (TIMES(I) + TIMES(I + 1)) / 2.
          IF(TMID .GT. RTIME) THEN
            ISTEP = I
            GO TO 33
          END IF
   32   CONTINUE
   33   WRITE(NOUT,3020)RTIME,TIMES(ISTEP),ISTEP
        GO TO 5
      END IF
C
   35 CONTINUE
C
C Step
C
      IF (KVALUE(2) .NE. 1 .AND. KVALUE(2) .NE. 2)THEN
        IF(CVALUE(2)(1:3) .EQ. 'ALL')THEN
          ISTEP = -1
          IDEF = 0
          WRITE(NOUT,3005)
          GO TO 40
        ELSE
          WRITE(NOUT,3015)CVALUE(2)
          GO TO 5
        END IF
      ELSE
        ISTEP = IVALUE(2)
        WRITE(NOUT,3025)ISTEP,TIMES(ISTEP)
        GO TO 5
      END IF
C
   39 CONTINUE
C
C output time
C
      IF (ISTEP .EQ. -1)THEN
        WRITE(NOUT,3040)
        GO TO 5
      END IF
      IF (KVALUE(2) .EQ. 0)THEN
        IF (KVALUE(3) .EQ. 1 .OR. KVALUE(3) .EQ. 2)THEN
          OUTTIM = RVALUE(3)
        END IF
      ELSE IF (KVALUE(2) .EQ. 1 .OR. KVALUE(2) .EQ. 2)THEN
        OUTTIM = RVALUE(2)
      ELSE
        WRITE(NOUT,3030)CVALUE(2),CVALUE(3)
      END IF
      GO TO 5
C
   40 CONTINUE
C
C List times
C
      WRITE (NOUT,4000)
      WRITE (NOUT,4010)(TIMES(I),I=1,NUMTIM)
      GO TO 5
C
   50 CONTINUE
C
C Scheme
C
      IF (KVALUE(2) .NE. 2)THEN
        WRITE(NOUT,5000)CVALUE(2)
        GO TO 5
      END IF
      ISCHEM = IVALUE(2)
      DO 55 I = 1, IMP
        MP(3,I) = ISCHEM
   55 CONTINUE
      IF (ISCHEM .EQ. 0)THEN
        WRITE(NOUT,5010)ISCHEM
      ELSE IF (ISCHEM .EQ. 1)THEN
        WRITE(NOUT,5020)ISCHEM
      ELSE IF (ISCHEM .EQ. 2)THEN
        WRITE(NOUT,5030)ISCHEM
      ELSE IF (ISCHEM .EQ. 3)THEN
        WRITE(NOUT,5035)ISCHEM
      ELSE
        WRITE(NOUT,5040)ISCHEM
      END IF
      GO TO 5
C
   60 CONTINUE
C
C Searchbox (tolerance)
C
      IF (KVALUE(2) .EQ. 1 .OR. KVALUE(2) .EQ. 2) THEN
        TOLSHL = RVALUE(2)
        TOLQAD = RVALUE(2)
        TOLHEX = RVALUE(2)
        WRITE(NOUT,6000)TOLSHL
      ELSE
        WRITE(NOUT,6010)CVALUE(2)
      END IF
      GO TO 5
C
   70 CONTINUE
C
C Deformed vs undeformed processing
C
      IF (KVALUE(2) .NE. 2)THEN
        WRITE(NOUT,7000)CVALUE(2)
        GO TO 5
      END IF
      IDEF = IVALUE(2)
      IF (IDEF .EQ. 0)THEN
        WRITE(NOUT,7010)idef
      ELSE IF (IDEF .EQ. 1 .OR. IDEF .EQ. 2)THEN
        IF (ISTEP .EQ. -1)THEN
          WRITE(NOUT,7020)idef
          IDEF = 0
        END IF
        IF (IDEF .EQ. 1)THEN
          WRITE(NOUT,7030)idef
        ELSE IF (IDEF .EQ. 2) THEN
          WRITE(NOUT,7035)idef
        END IF
      ELSE 
        WRITE(NOUT,7040)idef
      END IF
      GO TO 5
C
   80 CONTINUE
C
C Map definition - donor mesh e-block to recipient mesh e-block
C
      IF (KVALUE(2) .EQ. 0 .AND. CVALUE(2) .EQ. 'ALL')THEN
        IF (CVALUE(3) .NE. 'TO')THEN
          WRITE(NOUT,8010)
          WRITE(NTPOUT,8010)
          GO TO 5
        END IF
        IF (KVALUE(4) .NE. 1 .AND. KVALUE(4) .NE. 2)THEN
          WRITE(NOUT,8020)CVALUE(4)
          WRITE(NTPOUT,8020)CVALUE(4)
          GO TO 5
        END IF
        IF (KVALUE(5) .EQ. 0 .AND. CVALUE(5)(1:3) .EQ. 'SCH' .AND.
     &      (KVALUE(6) .EQ. 1 .OR. KVALUE(6) .EQ. 2))THEN
          IF (IVALUE(6) .GE. 0 .AND. IVALUE(6) .LE. 3) THEN
            ISCHEM = IVALUE(6)
          END IF
        END IF
C ... Check for valid id
        ID = IVALUE(4)
        idl = locint(id, nblksb, idb)
        if (idl .eq. 0) then
          write (nout, 8090) id
          go to 5
        end if
        DO 82 I = 1, NBLKSA
          MP(1,I) = IDA(I)
          MP(2,I) = IVALUE(4)
          MP(3,I) = ISCHEM
   82   CONTINUE
        IMP = NBLKSA
        GO TO 5
      ELSE IF (KVALUE(4) .EQ. 0 .AND. CVALUE(4) .EQ. 'ALL')THEN
        IF (CVALUE(3) .NE. 'TO')THEN
          WRITE(NOUT,8010)
          WRITE(NTPOUT,8010)
          GO TO 5
        END IF
        IF (KVALUE(2) .NE. 1 .AND. KVALUE(2) .NE. 2)THEN
          WRITE(NOUT,8030)CVALUE(2)
          WRITE(NTPOUT,8030)CVALUE(2)
          GO TO 5
        END IF
        IF (KVALUE(5) .EQ. 0 .AND. CVALUE(5)(1:3) .EQ. 'SCH' .AND.
     &      (KVALUE(6) .EQ. 1 .OR. KVALUE(6) .EQ. 2))THEN
          IF (IVALUE(6) .GE. 0 .AND. IVALUE(6) .LE. 3) THEN
            ISCHEM = IVALUE(6)
          END IF
        END IF

C ... Check for valid id
        ID = IVALUE(2)
        idl = locint(id, nblksa, ida)
        if (idl .eq. 0) then
          write (nout, 8090) id
          go to 5
        end if

        DO 84 I = 1, NBLKSB
          MP(2,I) = IDB(I)
          MP(1,I) = IVALUE(2)
          MP(3,I) = ISCHEM
   84   CONTINUE
        IMP = NBLKSB
        GO TO 5
      ELSE IF (KVALUE(2) .EQ. 1 .OR. KVALUE(2) .EQ. 2 .AND.
     &         KVALUE(4) .EQ. 1 .OR. KVALUE(4) .EQ. 2)THEN
        IF (CVALUE(3) .NE. 'TO')THEN
          WRITE(NOUT,8010)
          WRITE(NTPOUT,8010)
          GO TO 5
        END IF
        IF (KVALUE(5) .EQ. 0 .AND. CVALUE(5)(1:3) .EQ. 'SCH' .AND.
     &      (KVALUE(6) .EQ. 1 .OR. KVALUE(6) .EQ. 2))THEN
          IF (IVALUE(6) .GE. 0 .AND. IVALUE(6) .LE. 3) THEN
            ISCHEM = IVALUE(6)
          END IF
        END IF

C ... Check for valid ids
        ID = IVALUE(2)
        idl = locint(id, nblksa, ida)
        if (idl .eq. 0) then
          write (nout, 8090) id
          go to 5
        end if

        ID = IVALUE(4)
        idl = locint(id, nblksb, idb)
        if (idl .eq. 0) then
          write (nout, 8090) id
          go to 5
        end if

        IMP2 = IMP2 + 1
        IMP = IMP2
        MP(1,IMP2) = IVALUE(2)
        MP(2,IMP2) = IVALUE(4)
        MP(3,IMP2) = ISCHEM
        WRITE(NOUT,8040)
        DO 86 I = 1, IMP
          WRITE(NOUT,8050)MP(1,I),MP(2,I),MP(3,I)
   86   CONTINUE
        WRITE(NOUT,8070)
        DO 88 I = 1, NBLKSB
          WRITE(NOUT,8080)IDB(I)
   88   CONTINUE
        GO TO 5
      ELSE
        WRITE(NOUT,8060)CVALUE(2)
        WRITE(NTPOUT,8060)CVALUE(2)
        GO TO 5
      END IF
C
   89 CONTINUE
C
C Read integer flag for accuracy checks (comparison of
C various quantities between donor and recipient meshes
C
      IF (KVALUE(2) .NE. 2) THEN
        WRITE (NOUT,8900)CVALUE(2)
        WRITE (NTPOUT,8900)CVALUE(2)
        GO TO 5
      END IF
      IACCU = IVALUE(2)
      IF (IACCU .EQ. 0)THEN
        WRITE(NOUT,8910)IACCU
        WRITE(NTPOUT,8910)IACCU
      ELSE IF (IACCU .EQ. 1)THEN
        WRITE(NOUT,8920)IACCU
        WRITE(NTPOUT,8920)IACCU
      ELSE
        WRITE(NOUT,8930)IACCU
        WRITE(NTPOUT,8930)IACCU
      END IF
      GO TO 5
C
   90 CONTINUE
C
C Stop execution
C
      WRITE(NOUT,9000)
      WRITE(NTPOUT,9000)
      CALL ERROR('RDINPT','YOU ELECTED TO TERMINATE THE PROGRAM',' ',
     &           0,' ',0,' ',' ',1) 
C
  100 CONTINUE
C
C Continue execution (run)
C 
C sort map array (MP) on second entry (recipient mesh element block)
C this is required because of way mapping of multiple donor mesh
C element blocks into one recipient mesh element blocks is 
C implemented (required to have all such maps located sequentially)
C a simple sort on the 2nd entry accomplishes this and is easier
C than rewriting the offending algorythm
C
C
      IBOTOM = IMP - 1
 110  ISWICH = 1
      DO 120 I = 1, IBOTOM
        IF (MP(2,I) .LE. MP(2,I+1))THEN
          GO TO 120
        ELSE
          ITEMP = MP(1,I)
          MP(1,I) = MP(1,I+1)
          MP(1,I+1) = ITEMP
          ITEMP = MP(2,I)
          MP(2,I) = MP(2,I+1)
          MP(2,I+1) = ITEMP
          ITEMP = MP(3,I)
          MP(3,I) = MP(3,I+1)
          MP(3,I+1) = ITEMP
          ISWICH = I
          GO TO 120
        END IF
  120 CONTINUE
      IF (ISWICH .EQ. 1)THEN
        GO TO 130
      ELSE
        IBOTOM = ISWICH - 1
        GO TO 110
      END IF
  130 CONTINUE
C
C end sort
C
      WRITE(NOUT,10000)
      WRITE(NTPOUT,10000)
      IF (ISTEP .EQ. -1) THEN
        WRITE(NOUT,10010)
        WRITE(NTPOUT,10010)
      ELSE
        WRITE(NOUT,10020)RTIME,TIMES(ISTEP),ISTEP
        WRITE(NTPOUT,10020)RTIME,TIMES(ISTEP),ISTEP
      END IF
      IF (IDEF .EQ. 0)THEN
        WRITE(NOUT,7010)idef
        WRITE(NTPOUT,7010)idef
      ELSE IF (IDEF .EQ. 1)THEN
        WRITE(NOUT,7030)idef
        WRITE(NTPOUT,7030)idef
      ELSE IF (IDEF .EQ. 2) THEN
        WRITE(NOUT,7035)idef
        WRITE(NTPOUT,7035)idef
      END IF
      IF (ISCHEM .EQ. 0)THEN
         WRITE(NOUT,5010)ISCHEM
         WRITE(NTPOUT,5010)ISCHEM
      ELSE IF (ISCHEM .EQ. 2)THEN
         WRITE(NOUT,5030)ISCHEM
         WRITE(NTPOUT,5030)ISCHEM
      ELSE IF (ISCHEM .EQ. 3)THEN
         WRITE(NOUT,5035)ISCHEM
         WRITE(NTPOUT,5035)ISCHEM
      ELSE
         WRITE(NOUT,5020)ISCHEM
         WRITE(NTPOUT,5020)ISCHEM
      END IF
      WRITE(NOUT,10030)
      WRITE(NTPOUT,10030)
      DO 190 I = 1, IMP
        WRITE(NOUT,8050)MP(1,I),MP(2,I),MP(3,I)
        WRITE(NTPOUT,8050)MP(1,I),MP(2,I),MP(3,I)
  190 CONTINUE
      RETURN
 1000 FORMAT ('1',/////)
 1010 FORMAT (//,5X,'***************************************************
     1********',//)
 1020 FORMAT (5X,'VERSION -- ',A32,//)
 1030 FORMAT (5X,'DATE OF EXECUTION -- ',A32,//)
 1040 FORMAT (5X,'TIME OF EXECUTION -- ',A32,//)
 1050 FORMAT (2X,'READING INPUT DATA' ,//)
 1060 FORMAT (5X,'MAPVAR INPUT - SYNTAX:',/,' KEY_WORD <value>',//,
     1' HELp                  - REPEATS THIS MESSAGE',/,
     2' HELp <key_word>       - DESCRIPTION OF COMMAND - key_word',/,
     3' SCHeme <int>          - MAPPING SCHEME TO USE',/,
     4' DEFormed <int>        - ORIGINAL OR DEFORMED GEOMETRY',/,
     5' LISt TIMes            - WRITES A LIST OF TIMES AVAILABLE',/,
     6' TIMes <real or ALL>   - TIME TO BE MAPPED',/,
     7' STEp <int or ALL>     - TIME STEP TO BE MAPPED',/,
     8' OUTput TIMe <real>    - TIME TO BE WRITTEN TO OUTPUT FILE',/,
     9' SEArchbox <real>      - SEARCH BOX TO BE USED')
 1061 FORMAT (
     1' MAP <int or ALL> TO <int or ALL> SCHeme <int>',/,
     2'                       - ELEMENT BLOCK MAPPING ',/,
     3' CHEck <int>           - CHECK ACCURACY OF MAPPING',/,
     3' RUN                   - END INPUT - RESUME PROGRAM',/,
     4' STOP                  - TERMINATES THE PROGRAM')  
C
 1100 FORMAT(5X,'UNKNOWN INPUT - READING',A20,/,
     1'          PLEASE TRY AGAIN')
C
 2000 FORMAT(5X,'TIMe <real or ALL>',//,
     1'          IF A REAL NUMBER VALUE IS ENTERED, IT REPRESENTS',/,
     2'          THE TIME (STEP) SELECTED AT WHICH VARIABLES WILL',/,
     3'          BE MAPPED FROM THE DONOR TO THE RECIPIENT MESH',/,
     4'          DEFAULT - the last time step in restart file',/,
     5'          IF *ALL* IS ENTERED, ALL THE TIME STEPS IN THE',/,
     6'          DONOR MESH WILL BE MAPPED.',/,
     7'          NOTE: ONLY ORIGINAL GEOMETRY MAPPING IS ALLOWED',/,
     8'          *DEFORMED 0*')
 2010 FORMAT(5x,'LISt TIMes',//,
     1'          LIST TIMES COMMAND READS THE DONOR MESH FILE',/,
     2'          AND ECHOS TO THE CRT A LIST OF TIMES AVAILABLE TO',/,
     3'          THE USER.',//,
     4'          DEFAULT - none')
 2020 FORMAT(5x,'SCHeme <int>',//,
     1'          SCHEME COMMAND PICKS THE MAPPING SCHEME TO USE',/,
     2'          SCHEME 0 - NODAL BASED SIMPLE INTERPOLATION',/,
     3'          SCHEME 1 - NODAL BASED LEAST SQUARES',/,
     4'          SCHEME 2 - DIRECT TRANSFER',/,
     5'          SCHEME 3 - ELEMENT CENTROID BASED LEAST SQUARES',/,
     6'          DEFAULT  - 1')
 2030 FORMAT(5X,'DEFormed <int>',//,
     1'          THE DEFORMED COMMAND SELECTS USE OF ORIGINAL OR',/,
     2'          DEFORMED GEOMETRY',/,
     3'          DEFORMED 0 - ORIGINAL GEOMETRY - COORDINATES ARE',/,
     4'                       NOT MODIFIED BY DISPLACEMENTS',/,
     5'                       RECIPIENT MESH BOUNDARIES IDENTICAL',/,
     6'                       TO UNDEFORMED DONOR MESH BOUNDARIES')
 2040 FORMAT(10X,'DEFORMED 1 - DEFORMED GEOMETRY - DISPLACEMENTS ARE',/,
     1'                       ADDED TO DONOR MESH COORDINATES PRIOR',/,
     2'                       TO USE AND THE MAPPED DISPLACEMENTS',/,
     3'                       ARE SUBTRACTED FROM RECIPIENT MESH',/,
     4'                       COORDINATES PRIOR TO OUTPUT',/,
     5'                       RECIPIENT MESH BOUNDARIES IDENTICAL',/,
     6'                       TO DEFORMED DONOR MESH BOUNDARIES')
 2045 FORMAT(10X,'DEFORMED 2 - MESH ANNEALING - MAP IS PERFORMED IN',/,
     1'                       DEFORMED COORDINATES (DISPLACEMENTS',/,
     2'                       ARE ADDED TO DONOR MESH COORDINATES',/,
     3'                       PRIOR TO USE. OUTPUT DISPLACEMENTS',/,
     4'                       ARE SET TO ZERO.',/,
     5'                       THIS OPTION WAS REQUESTED FOR GOMA',/,
     6'          DEFAULT - 1')
 2050 FORMAT(10X,'SEArchbox <REAL>',//,
     1'           THE SEARCHBOX IS A BOUNDING BOX AROUND THE DONOR',/,
     2'           MESH ELEMENT OR GROUP OF ELEMENTS IN WHICH THE',/,
     3'           SEARCH ROUTINE ATTEMPTS TO FIND A POINT',/,
     4'           (RECIPIENT MESH NODE OR ELEMENT CENTROID).',/,
     5'           THE SEARCH BOX CAN BE INCREASED MODESTLY',/,
     5'           TO OVERCOME MESHING ISSUES - AS DIFFERENT',/,
     6'           DISCRETIZATION ON A RADIUS.',/,
     7'           BE GENTLE - INCREASE CAN HAVE A LARGE EFFECT ON',/,
     8'           RUN TIMES. GREATER THAN 1.0 IS NOT RECOMMENDED.',/,
     9'           DEFAULT = 0.01')
 2060 FORMAT(10X,'MAP <int or ALL> TO <int or ALL> SCHEME <INT>',//,
     1'           THE MAP COMMAND ALLOWS THE USER TO DEFINE THE',/,
     2'           DONOR MESH ELEMENT BLOCK I.D. TO BE MAPPED INTO',/,
     3'           THE RECIPIENT MESH ELEMENT BLOCK I.D. ALONG WITH',/, 
     4'           THE SCHEME TO USE FOR THIS MAPPING. IF THE',/,
     5'           VALUE "ALL" IS ENTERED THE OTHER FIELD MUST',/,
     6'           CONTAIN AN INTEGER',/,
     7'           DEFAULT:',/,
     8'           1 TO 1   2 TO 2   ETC',/,
     9'           DEFAULT SCHEME value eneterd with SCH command or 1')
 2070 FORMAT(10X,'CHEck <int>',//,
     1'           THE CHECK COMMAND COMPUTES QUANTITIES FOR',/,
     2'           COMPARISON BETWEEN THE DONOR AND RECIPIENT',/,
     3'           MESHES',/,
     4'           0 - NO CHECK QUANTITIES COMPUTED',/,
     5'           1 - ALL APPROPRIATE QUANTITIES COMPUTED',/,
     6'           DEFAULT - 0')
C
 3000 FORMAT(5X,'TIME YOU HAVE ENTERED - TIMES ALL',/,
     1'          ALL THE TIME STEPS WILL BE MAPPED',/,
     2'          ONLY UNDEFORMED GEOMETRY PROCESSING IS',/,
     3'          IMPLEMENTED WITH *TIMES ALL* TIME STEP INPUT',/,
     4'          DEFORMED 0',///)
 3005 FORMAT(5X,'STEP YOU HAVE ENTERED - STEPS ALL',/,
     1'          ALL THE TIME STEPS WILL BE MAPPED',/,
     2'          ONLY UNDEFORMED GEOMETRY PROCESSING IS',/,
     3'          IMPLEMENTED WITH *STEPS ALL* TIME STEP INPUT',/,
     4'          DEFORMED 0',///)
 3010 FORMAT(5X,'READING TIMES COMMAND',/,
     1'          EXPECTED A REAL NUMBER IN FIELD 2',/,
     2'          READ',A20,/,
     3'          PLEASE TRY AGAIN')
 3015 FORMAT(5X,'READING STEP COMMAND',/,
     1'          EXPECTED A REAL NUMBER IN FIELD 2',/,
     2'          READ',A20,/,
     3'          PLEASE TRY AGAIN')
 3020 FORMAT(5X,'TIME YOU HAVE ENTERED',F14.6,/,
     1'          CLOSEST TIME ON DATABASE',F14.6,/,
     2'          TIME STEP',I5,///)
 3025 FORMAT(5X,'STEP YOU HAVE ENTERED',I3,/,
     1'          TIME ',F14.6,///)
 3030 FORMAT(5X,'READING OUTPUT TIME COMMAND',/,
     1'          EXPECTED A REAL NUMBER IN FIELD 2 OR 3',/,
     2'          READ',A20,A20,/,
     3'          PLEASE TRY AGAIN')
 3040 FORMAT(5X,'YOU HAVE ALREADY SELECTED TO PROCESS',/,
     1'          ALL TIME STEPS AVAILABLE. YOU CANNOT ALSO',/,
     2'          CHANGE THE OUTPUT TIME')
C
 4000 FORMAT(5X,'TIMES AVAILABLE FROM THE RESTART FILE')
 4010 FORMAT(5X,/,E14.6)
C
 5000 FORMAT(5X,'READING SCHEME COMMAND',/,
     1'          EXPECTED AN INTEGER IN FIELD 2',/,
     2'          READ',A20,/,
     3'          PLEASE TRY AGAIN')
 5010 FORMAT(//,5X,'YOU HAVE ENTERED SCHEME ',I5,/,
     1'          NODAL BASED SIMPLE INTERPOLATION')
 5020 FORMAT(//,5X,'YOU HAVE ENTERED SCHEME ',I5,/,
     1'          NODAL BASED LEAST SQUARES')
 5030 FORMAT(//,5X,'YOU HAVE ENTERED SCHEME ',I5,/,
     1'          DIRECT TRANSFER')
 5035 FORMAT(//,5X,'YOU HAVE ENTERED SCHEME ',I5,/,
     1'          ELEMENT CENTROID BASED LEAST SQUARES')
 5040 FORMAT(5X,'YOU HAVE ENTERED SCHEME ',I5,/,
     1'          THIS SCHEME HAS NOT BEEN IMPLEMENTED',/,
     2'          PLEASE TRY AGAIN')
C
 6000 FORMAT(//,5X,'YOU HAVE ENTERED SEARCH ',F12.4,/,
     1'          VALUES GREATER THAN 1. ARE NOT RECOMMENDED')
 6010 FORMAT(5X,'READING SEArch COMMAND',/,
     1'          EXPECTED A REAL NUMBER IN FIELD 2',/,
     2'          READ',A20,/,
     3'          PLEASE TRY AGAIN')
C
 7000 FORMAT(5X,'READING DEFORMED GEOMETRY COMMAND',/,
     1'          EXPECTED AN INTEGER IN FIELD 2',/,
     2'          READ',A20,/,
     3'          PLEASE TRY AGAIN')
 7010 FORMAT(//,5X,'YOU HAVE ENTERED - DEFORMED',I5,/,
     1'          ORIGINAL GEOMETRY - COORDINATES ARE NOT',/,
     2'          MODIFIED BY DISPLACEMENTS')
 7020 FORMAT(5X,'YOU HAVE ENTERED - DEFORMED',I5,/,
     1'          ONLY ORIGINAL GEOMETRY PROCESSING IS',/,
     2'          COMPATIBLE WITH *TIMES ALL* TIME STEP INPUT',/,
     3'          DEFORMED 0',///)
 7030 FORMAT(//,5X,'YOU HAVE ENTERED - DEFORMED',I5,/,
     1'          DEFORMED GEOMETRY - DISPLACEMENTS ARE',/,
     2'          ADDED TO DONOR MESH COORDINATES AND',/,
     3'          SUBTRACTED FROM RECIPIENT MESH COORDINATES')
 7035 FORMAT(//,5X,'YOU HAVE ENTERED - DEFORMED',I5,/,
     1'          MESH ANNEALING - DISPLACEMENTS ARE',/,
     2'          ADDED TO DONOR MESH COORDINATES.',/,
     3'          DISPLACEMENTS ARE ZERO ON THE INTERPOLATED MESH.')
 7040 FORMAT(//,5X,'YOU HAVE ENTERED - DEFORMED',I5,/,
     1'          value - MUST BE EITHER 0 OR 1',/,
     2'          PLEASE TRY AGAIN')
C
 8010 FORMAT(5X,'READING MAP COMMAND',/,
     &'          EXPECTING THE CHARACTER STRING "TO" ',/,
     &'          COMMAND SYNTAX IS:',/,
     &'          MAP <int or "ALL"> TO <int or "ALL">',//,
     &'          PLEASE TRY AGAIN')
 8020 FORMAT(5X,'READING MAP COMMAND',/,
     &'          EXPECTING AN INTEGER AFTER "MAP ALL TO"',/,
     &'          READ',A20,/,
     &'          PLEASE TRY AGAIN')
 8030 FORMAT(5X,'READING MAP COMMAND',/,
     &'          EXPECTING AN INTEGER IN FIELD 2 OF COMMAND',/,
     &'          "MAP ?? TO ALL"',/,
     &'          READ',A20,/,
     &'          PLEASE TRY AGAIN')
 8040 FORMAT(5X,'READING MAP COMMAND',/,
     &'          MAP AS ENTERED SO FAR IS:')
 8050 FORMAT(5X,'MAP ',I5,' TO', I5, ' SCHEME', I5)
 8060 FORMAT(5X,'READING MAP COMMAND',/,
     &'          EXPECTING EITHER AN INTEGER OR THE STRING "ALL"',/,
     &'          READ',A20,/,
     &'          PLEASE TRY AGAIN')
 8070 FORMAT(//,5X,'AVAILABLE RECIPIENT MESH ELEMENT BLOCK I.D.')
 8080 FORMAT(10X,I7)
 8090 FORMAT(//,5x,'ERROR: The entered id ', i5,
     *  ' is not a valid block id.',/)
C
 8900 FORMAT(5X,'READING CHECK ACCURACY COMMAND',/,
     1'          EXPECTED AN INTEGER IN FIELD 2',/,
     2'          READ',A20,/,
     3'          PLEASE TRY AGAIN')
 8910 FORMAT(//,5X,'YOU HAVE ENTERED - CHECK',I5,/,
     1'          NO ACCURACY CHECKING WILL BE DONE')
 8920 FORMAT(//,5X,'YOU HAVE ENTERED - CHECK',I5,/,
     1'          ALL APPROPRIATE QUANTITIES WILL BE COMPUTED',/,
     2'          AND WRITTEN TO THE TEXT OUTPUT FILE FOR',/,
     3'          COMPARISON BETWEEN THE DONOR AND RECIPIENT MESHES')
 8930 FORMAT(//,5X,'YOU HAVE ENTERED - CHECK',I5,/,
     1'          ONLY VALUES 0 OR 1 HAVE BEEN IMPLEMENTED')
C
 9000 FORMAT(5X,'YOU HAVE ELECTED TO TERMINATE THE PROGRAM',/,
     1'          NOTHING WILL BE COMPUTED OR SAVED')
C
10000 FORMAT(5X,'LEAVING RDINPT - VALUES USED ARE:',//)
10010 FORMAT(//5X,'YOU HAVE ENTERED FOR THE TIMES COMMAND',/,
     1'          *TIMES ALL* - ALL TIME STEPS WILL BE MAPPED',//)
10020 FORMAT(//5X,'YOU HAVE ENTERED FOR THE TIMES COMMAND',/,
     1'          rtime = ',F14.6,/,
     2'          CLOSEST TIME ON DATABASE',/,
     3'          ctime = ',f14.6,/,
     4'          TIME STEP',/,
     5'          istep =',i5,//)
10030 FORMAT(5X,'MAP TO BE USED:')
C
      END

C=======================================================================
      INTEGER FUNCTION LOCINT (INT, LENLST, INTLST)
C=======================================================================
C   --*** LOCINT *** (ETCLIB) Find integer in list
C   --   Written by Amy Gilkey - revised 11/10/87
C   --
C   --LOCINT returns the index of the given integer in a list of integers.
C   --If the integer is not in the list, returns 0.
C   --
C   --Parameters:
C   --   INT - IN - the integer to be searched for
C   --   LENLST - IN - the number of integers in the list
C   --   INTLST - IN - the list of integers to be searched

      INTEGER INT
      INTEGER LENLST
      INTEGER INTLST(*)

      DO 10 LOCINT = 1, LENLST
         IF (INT .EQ. INTLST(LOCINT)) GOTO 20
   10 CONTINUE
      LOCINT = 0

   20 CONTINUE
      RETURN
      END





