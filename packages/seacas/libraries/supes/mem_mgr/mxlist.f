C    Copyright(C) 2008 Sandia Corporation.  Under the terms of Contract
C    DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C    certain rights in this software
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C    * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C              
C    * Redistributions in binary form must reproduce the above
C      copyright notice, this list of conditions and the following
C      disclaimer in the documentation and/or other materials provided
C      with the distribution.
C                            
C    * Neither the name of Sandia Corporation nor the names of its
C      contributors may be used to endorse or promote products derived
C      from this software without specific prior written permission.
C                                                    
C    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C    A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C    
      SUBROUTINE MXLIST (UNIT, OFFSET,
     *   DICT, DPOINT, LDICT, NNAMES,
     *   VOID, LVOID, NVOIDS, CHRCOL)
C
      IMPLICIT INTEGER (A-Z)
C
C     This routine lists the contents of the tables of the
C     memory manager.
C
C***********************************************************************
C
C     UNIT     Output unit number
C     OFFSET   Offset to internal reference vector
C     DICT     Dictionary table
C     DPOINT   Pointer table
C     LDICT    Dimension of dictionary
C     NNAMES   Number of names in the dictionary
               CHARACTER*8 DICT(LDICT,CHRCOL)
               DIMENSION DPOINT(LDICT,CHRCOL,3), NNAMES(2)
C     VOID     Void table
C     LVOID    Dimension of void table
C     NVOIDS   Number of voids
               DIMENSION VOID(LVOID,CHRCOL,2), NVOIDS(2)
C     CHRCOL   Column number for character tables
C
C***********************************************************************
C
      CHARACTER*8 TNAME1, TNAME2
C
      TOFF = OFFSET
      DO 290 ICOL = 1, CHRCOL
         IF (ICOL .EQ. 2) THEN
            TOFF = 0
            WRITE (UNIT, 10090)
         END IF
C        DICTIONARY.
C
         MNDICT = 0
         MCDICT = 0
         DO 100 I = 1, NNAMES(ICOL)
            MNDICT = MNDICT + MAX(0,DPOINT(I,ICOL,2))
            MCDICT = MCDICT + MAX(0,DPOINT(I,ICOL,3))
  100    CONTINUE
         WRITE (UNIT,10000)
         DO 110 I = 1, NNAMES(ICOL)
            IF (DPOINT(I,ICOL,2) .GE. 0) THEN
               WRITE (UNIT,10010) I, DICT(I,ICOL),
     *            DPOINT(I,ICOL,1)+TOFF,
     *            DPOINT(I,ICOL,2), DPOINT(I,ICOL,3)
            ELSE
               WRITE (UNIT,10020) I, DICT(I,ICOL),
     *            DPOINT(I,ICOL,2), DPOINT(I,ICOL,3)
            END IF
  110    CONTINUE
         WRITE (UNIT,10040) MNDICT, MCDICT
C
C        VOID TABLE.
C
         MNVOID = 0
         DO 120 I = 1, NVOIDS(ICOL)
            MNVOID = MNVOID + VOID(I,ICOL,2)
  120    CONTINUE
         WRITE (UNIT,10070)
         WRITE (UNIT,10060) (I,VOID(I,ICOL,1)+TOFF,
     *      VOID(I,ICOL,2),I=1,NVOIDS(ICOL))
         WRITE (UNIT,10050) MNVOID
C
C        OUTPUT ORDERED LIST OF TABLES.
C
C        First sort dictionary into location order with unresolved
C        allocations first.
C
         JSTRT = 2
         DO 150 I = 1, NNAMES(ICOL)-1
            IF (DPOINT(I,ICOL,2) .GE. 0) THEN
               DO 130 J = MAX(JSTRT,I+1), NNAMES(ICOL)
                  IF (DPOINT(J,ICOL,2) .LT. 0) THEN
                     JSTRT = J+1
                     TNAME1 = DICT(I,ICOL)
                     TNAME2 = DICT(J,ICOL)
                     DICT(I,ICOL) = TNAME2
                     DICT(J,ICOL) = TNAME1
                     TEMP = DPOINT(I,ICOL,1)
                     DPOINT(I,ICOL,1) = DPOINT(J,ICOL,1)
                     DPOINT(J,ICOL,1) = TEMP
                     TEMP = DPOINT(I,ICOL,2)
                     DPOINT(I,ICOL,2) = DPOINT(J,ICOL,2)
                     DPOINT(J,ICOL,2) = TEMP
                     TEMP = DPOINT(I,ICOL,3)
                     DPOINT(I,ICOL,3) = DPOINT(J,ICOL,3)
                     DPOINT(J,ICOL,3) = TEMP
                     GO TO 140
                  END IF
  130          CONTINUE
               GO TO 160
            END IF
  140       CONTINUE
  150    CONTINUE
  160    CONTINUE
         DO 170 NDEFER = 1, NNAMES(ICOL)
            IF (DPOINT(NDEFER,ICOL,2) .GE. 0) GO TO 180
  170    CONTINUE
  180    NDEFER = NDEFER - 1
         DO 200 I = NDEFER+1, NNAMES(ICOL)-1
            DO 190 J = I+1, NNAMES(ICOL)
               IF (DPOINT(J,ICOL,1) .LT. DPOINT(I,ICOL,1)) THEN
                  TNAME1 = DICT(I,ICOL)
                  TNAME2 = DICT(J,ICOL)
                  DICT(I,ICOL) = TNAME2
                  DICT(J,ICOL) = TNAME1
                  TEMP = DPOINT(I,ICOL,1)
                  DPOINT(I,ICOL,1) = DPOINT(J,ICOL,1)
                  DPOINT(J,ICOL,1) = TEMP
                  TEMP = DPOINT(I,ICOL,2)
                  DPOINT(I,ICOL,2) = DPOINT(J,ICOL,2)
                  DPOINT(J,ICOL,2) = TEMP
                  TEMP = DPOINT(I,ICOL,3)
                  DPOINT(I,ICOL,3) = DPOINT(J,ICOL,3)
                  DPOINT(J,ICOL,3) = TEMP
               END IF
  190       CONTINUE
  200    CONTINUE
C
C        STARTING STUFF FOR LOOP
C
         DO 210 IDICT = 1, NNAMES(ICOL)
            IF (DPOINT(IDICT,ICOL,2) .GE. 0) GO TO 220
  210    CONTINUE
  220    CONTINUE
         IF (IDICT .LE. NNAMES(ICOL) .AND. NVOIDS(ICOL) .GT. 0) THEN
            NXTLOC = MIN(DPOINT(IDICT,ICOL,1), VOID(1,ICOL,1))
         ELSE IF (IDICT .LE. NNAMES(ICOL)) THEN
            NXTLOC = DPOINT(IDICT,ICOL,1)
         ELSE IF (NVOIDS(ICOL) .GT. 0) THEN
            NXTLOC = VOID(1,ICOL,1)
         END IF
         IVOID = 1
C
         WRITE (UNIT, 10080)
         ILIST = 0
C
C        Deferred space names first.
C
         DTOT = 0
         DCTOT = 0
         DO 230 IDICT = 1, NNAMES(ICOL)
            IF (DPOINT(IDICT,ICOL,2) .LT. 0) THEN
               ILIST = ILIST + 1
               WRITE (UNIT,10020) ILIST, DICT(IDICT,ICOL),
     *            DPOINT(IDICT,ICOL,2),
     *            DPOINT(IDICT,ICOL,3)
               DTOT = DTOT - DPOINT(IDICT,ICOL,2)
               DCTOT = DCTOT + MAX (0, DPOINT(IDICT,ICOL,3))
            ELSE
               GO TO 240
            END IF
  230    CONTINUE
  240    CONTINUE
         IF (IDICT .GT. 1) WRITE (UNIT, 10030)
     *      'DEFERRED TOTAL', DTOT, DCTOT
C
C        LOOP
C
         TOTAL = 0
         SUBTOT = 0
         CSTOT = 0
         CTOT = 0
  250    CONTINUE
            IF (IDICT .LE. NNAMES(ICOL)
     *         .AND. IVOID .LE. NVOIDS(ICOL)) THEN
               IF (NXTLOC .EQ. DPOINT(IDICT,ICOL,1)) THEN
                  ILIST = ILIST + 1
                  WRITE (UNIT,10010) ILIST, DICT(IDICT,ICOL),
     *               DPOINT(IDICT,ICOL,1)+TOFF, DPOINT(IDICT,ICOL,2),
     *               DPOINT(IDICT,ICOL,3)
                  NXTLOC = NXTLOC + DPOINT(IDICT,ICOL,2)
                  SUBTOT = SUBTOT + DPOINT(IDICT,ICOL,2)
                  TOTAL = TOTAL + DPOINT(IDICT,ICOL,2)
                  CSTOT = CSTOT + MAX(0,DPOINT(IDICT,ICOL,3))
                  CTOT = CTOT + MAX(0,DPOINT(IDICT,ICOL,3))
                  IDICT = IDICT + 1
               ELSE IF (NXTLOC .EQ. VOID(IVOID,ICOL,1)) THEN
                  ILIST = ILIST + 1
                  WRITE (UNIT,10010) ILIST, ' ',
     *               VOID(IVOID,ICOL,1)+TOFF,
     *               VOID(IVOID,ICOL,2)
                  NXTLOC = NXTLOC + VOID(IVOID,ICOL,2)
                  SUBTOT = SUBTOT + VOID(IVOID,ICOL,2)
                  TOTAL = TOTAL + VOID(IVOID,ICOL,2)
                  IVOID = IVOID + 1
               ELSE
                  NXTLOC = MIN( DPOINT(IDICT,ICOL,1),
     *               VOID(IVOID,ICOL,1) )
                  WRITE (UNIT, 10030) 'BLOCK SIZE', SUBTOT, CSTOT
                  SUBTOT = 0
               END IF
            ELSE IF (IDICT .LE. NNAMES(ICOL)) THEN
               IF (NXTLOC .EQ. DPOINT(IDICT,ICOL,1)) THEN
                  ILIST = ILIST + 1
                  WRITE (UNIT,10010) ILIST, DICT(IDICT,ICOL),
     *               DPOINT(IDICT,ICOL,1)+TOFF, DPOINT(IDICT,ICOL,2),
     *               DPOINT(IDICT,ICOL,3)
                  NXTLOC = NXTLOC + DPOINT(IDICT,ICOL,2)
                  SUBTOT = SUBTOT + DPOINT(IDICT,ICOL,2)
                  TOTAL = TOTAL + DPOINT(IDICT,ICOL,2)
                  CSTOT = CSTOT + MAX(0,DPOINT(IDICT,ICOL,3))
                  CTOT = CTOT + MAX(0,DPOINT(IDICT,ICOL,3))
                  IDICT = IDICT + 1
               ELSE
                  NXTLOC = DPOINT(IDICT,ICOL,1)
                  WRITE (UNIT, 10030) 'BLOCK SIZE', SUBTOT, CSTOT
                  SUBTOT = 0
               END IF
            ELSE IF (IVOID .LE. NVOIDS(ICOL)) THEN
               IF (NXTLOC .EQ. VOID(IVOID,ICOL,1)) THEN
                  ILIST = ILIST + 1
                  WRITE (UNIT,10010) ILIST, ' ',
     *               VOID(IVOID,ICOL,1)+TOFF,
     *               VOID(IVOID,ICOL,2)
                  NXTLOC = NXTLOC + VOID(IVOID,ICOL,2)
                  SUBTOT = SUBTOT + VOID(IVOID,ICOL,2)
                  TOTAL = TOTAL + VOID(IVOID,ICOL,2)
                  IVOID = IVOID + 1
               ELSE
                  NXTLOC = VOID(IVOID,ICOL,1)
                  WRITE (UNIT, 10030) 'BLOCK SIZE', SUBTOT, CSTOT
                  SUBTOT = 0
               END IF
            ELSE
               GO TO 260
            END IF
         GO TO 250
  260    CONTINUE
         WRITE (UNIT, 10030) 'BLOCK SIZE', SUBTOT, CSTOT
         WRITE (UNIT, 10030) 'ALLOCATED TOTAL', TOTAL, CTOT
         WRITE (UNIT, 10030) '    GRAND TOTAL', TOTAL+DTOT, CTOT+DCTOT
C
C        SORT DICTIONARY BACK INTO NAME ORDER
C
         DO 280 I = 1, NNAMES(ICOL)-1
            DO 270 J = I+1, NNAMES(ICOL)
               IF (DICT(J,ICOL) .LT. DICT(I,ICOL)) THEN
                  TNAME1 = DICT(I,ICOL)
                  TNAME2 = DICT(J,ICOL)
                  DICT(I,ICOL) = TNAME2
                  DICT(J,ICOL) = TNAME1
                  TEMP = DPOINT(I,ICOL,1)
                  DPOINT(I,ICOL,1) = DPOINT(J,ICOL,1)
                  DPOINT(J,ICOL,1) = TEMP
                  TEMP = DPOINT(I,ICOL,2)
                  DPOINT(I,ICOL,2) = DPOINT(J,ICOL,2)
                  DPOINT(J,ICOL,2) = TEMP
                  TEMP = DPOINT(I,ICOL,3)
                  DPOINT(I,ICOL,3) = DPOINT(J,ICOL,3)
                  DPOINT(J,ICOL,3) = TEMP
               END IF
  270       CONTINUE
  280    CONTINUE
  290 CONTINUE
      RETURN
10000 FORMAT(//1X,50('*')//,7('* '),
     *   ' D I C T I O N A R Y  ',7('* ')//
     *   T33, 'NUMERIC', T43, 'CHARACTER'/
     *   T10,'NAME      LOCATION      LENGTH      LENGTH'/)
10010 FORMAT((1X,I4,2X,A8,3(2X,I10)))
10020 FORMAT((1X,I4,2X,A8,12X,2(2X,I10)))
10030 FORMAT (1X, A, T30, I10, 2X, I10/)
10040 FORMAT(/,T10,'TOTAL',T30,I10, 2X, I10)
10050 FORMAT(/,T10,'TOTAL',T24,I10)
10060 FORMAT((1X,I6,2(3X,I10)))
10070 FORMAT(///' * * * V O I D   T A B L E * * *'//
     *   '           LOCATION       LENGTH'/)
10080 FORMAT(//1X,50('*')/'0 ',6('* '),
     *   ' O R D E R E D   L I S T  ',6('* ')//
     *   ,T33, 'NUMERIC', T43, 'CHARACTER'/
     *   T10,'NAME      LOCATION      LENGTH      LENGTH'/)
10090 FORMAT(//1X,79('*')//,' SPLIT CHARACTER STORAGE'//1X,79('*')//)
      END
