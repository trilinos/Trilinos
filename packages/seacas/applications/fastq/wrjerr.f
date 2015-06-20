C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Government retains certain rights in this software.
C    
C    Redistribution and use in source and binary forms, with or without
C    modification, are permitted provided that the following conditions are
C    met:
C    
C        * Redistributions of source code must retain the above copyright
C          notice, this list of conditions and the following disclaimer.
C    
C        * Redistributions in binary form must reproduce the above
C          copyright notice, this list of conditions and the following
C          disclaimer in the documentation and/or other materials provided
C          with the distribution.
C    
C        * Neither the name of Sandia Corporation nor the names of its
C          contributors may be used to endorse or promote products derived
C          from this software without specific prior written permission.
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

C $Id: wrjerr.f,v 1.3 1998/07/14 18:20:17 gdsjaar Exp $
C $Log: wrjerr.f,v $
C Revision 1.3  1998/07/14 18:20:17  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.2  1991/04/10 19:57:37  gdsjaar
C Problem with internal read, kludged
C
c Revision 1.1.1.1  1990/11/30  11:17:55  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:17:54  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]WRJERR.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE WRJERR (MS, MR, NPNODE, NPELEM, MXNFLG, MXSFLG, NPREGN,
     &   NPNBC, NPSBC, IUNIT, NNN, KKK, NNXK, NODES, NELEMS, NNFLG,
     &   NNPTR, NNLEN, NSFLG, NSPTR, NSLEN, NVPTR, NVLEN, NSIDEN,
     &   MAPDXG, XN, YN, NXK, MAT, MAPGXD, MATMAP, NBCNOD, NNLIST,
     &   NBCSID, NSLIST, NVLIST, NUMMAT, LINKM, TITLE, ERR, EIGHT, NINE)
C***********************************************************************
C
C  SUBROUTINE WRJERR = WRITES JOE'S ERROR DATABASE MESH OUTPUT FILE
C
C***********************************************************************
C
      DIMENSION XN (NPNODE), YN (NPNODE), NXK (NNXK, NPELEM)
      DIMENSION MAT (NPELEM)
      DIMENSION NODES (NPNBC), NELEMS (NPSBC), NSIDEN (NPSBC)
      DIMENSION NNFLG (MXNFLG), NNLEN (MXNFLG), NNPTR (MXNFLG)
      DIMENSION NSFLG (MXSFLG), NSLEN (MXSFLG), NSPTR (MXSFLG)
      DIMENSION NVLEN (MXSFLG), NVPTR (MXSFLG), LINKM (2,  (MS+MR))
      DIMENSION MAPDXG (NPNODE), MAPGXD (NPNODE), MATMAP (3, NPREGN)
      DIMENSION IHOLD (9)
C
      CHARACTER*72 TITLE, DUMMY
C
      LOGICAL ERR, EIGHT, NINE, FOUND
C
      ERR = .TRUE.
C
C  WRITE OUT HEADER TITLE AND INFORMATION
C
      CALL INQSTR ('TITLE: ',TITLE)
      WRITE (IUNIT, 10000, ERR = 290)TITLE
      WRITE (IUNIT, 10010, ERR = 290)
C
C  WRITE OUT NODE BLOCK
C
      WRITE (IUNIT, 10020, ERR = 290)
      Z = 0.
      DO 100 I = 1, NNN
         WRITE (IUNIT, 10030, ERR = 290)I, XN (I), YN (I), Z
  100 CONTINUE
C
C  WRITE OUT ELEMENT BLOCKS
C
      DO 130 I = 1, NUMMAT
         CALL GETDUM (MATMAP (1, I), DUMMY, LEN)
         WRITE (IUNIT, 10040, ERR = 290)
         IF (NXK (3, MATMAP (2, I)) .EQ. 0) THEN
            INODE = 2
         ELSEIF (NXK (4, MATMAP (2, I)) .EQ. 0) THEN
            INODE = 3
         ELSEIF (EIGHT) THEN
            INODE = 8
         ELSEIF (NINE) THEN
            INODE = 9
         ELSE
            INODE = 4
         ENDIF
         IF (INODE .EQ. 8) THEN
            DO 110 KELEM = MATMAP (2, I), MATMAP (3, I)
               K = MAPGXD(KELEM)
               WRITE (IUNIT, 10050, ERR = 290) K, NXK (1, KELEM),
     &            NXK (3, KELEM), NXK (5, KELEM), NXK (7, KELEM),
     &            NXK (2, KELEM), NXK (4, KELEM), NXK (6, KELEM),
     &            NXK (8, KELEM)
  110       CONTINUE
         ELSE
            DO 120 KELEM = MATMAP (2, I), MATMAP (3, I)
               K = MAPGXD(KELEM)
               WRITE (IUNIT, 10060, ERR = 290) K, (NXK (J, KELEM),
     &            J = 1, INODE)
  120       CONTINUE
         ENDIF
  130 CONTINUE
C
C  WRITE OUT THE NODAL BOUNDARY CONDITIONS
C
      IF (NBCNOD.GT.0) THEN
         WRITE (IUNIT, 10070)
         DO 150 I = 1, NBCNOD
            J1 = NNPTR (I)
            J2 = NNPTR (I)+NNLEN (I)-1
            CALL GETDUM (NNFLG (I), DUMMY, LEN)
            WRITE (*, 10080) NNFLG(I)
            CALL INQSTR ('FIXED DEGREE OF FREEDOM: ',TITLE)
            READ (TITLE, '(I10)') IFIX
            DO 140 J = J1, J2
               WRITE (IUNIT, 10100, ERR = 290) NODES (J), IFIX
  140       CONTINUE
  150    CONTINUE
      ENDIF
C
C  WRITE OUT THE SIDE BOUNDARY FLAGS
C
      WRITE (IUNIT, 10120, ERR = 290)
      IF (NBCSID.GT.0) THEN
         WRITE (IUNIT, 10130, ERR = 290)
C         CALL MESAGE ('ELEMENT NUMBERING IS WRITTEN WITH ELEMENT' //
C     &      BOUNDARY FLAGS')
C         CALL INQTRU ('WOULD YOU LIKE TO CHANGE THIS TO NODES', IANS)
C         IF (IANS) THEN
C            DO 190 I = 1, NBCSID
C               J1 = NVPTR (I)
C               J2 = NVPTR (I) + NVLEN (I)-1
C               CALL GETDUM (NSFLG (I), DUMMY, LEN)
C               WRITE (IUNIT, 180, ERR = 200) DUMMY (1:LEN)
C               WRITE (IUNIT, 160, ERR = 200) (NSIDEN (J), J = J1, J2)
C 190        CONTINUE
C         ELSE
         DO 280 I = 1, NBCSID
            J1 = NSPTR (I)
            J2 = NSPTR (I)+NSLEN (I)-1
            CALL GETDUM (NSFLG (I), DUMMY, LEN)
            WRITE (*, 10090) NSFLG(I)
            CALL INQSTR ('PRESSURE MAGNITUDE: ',TITLE)
            READ (TITLE, '(F10.0)') PMAG
C
C  WRITE OUT THE SIDE 1 ELEMENTS
C
            FOUND = .FALSE.
            JHOLD = 0
            DO 170 J = J1, J2
               JJ1 = NSIDEN ( (J * 2) - 1)
               JJ2 = NSIDEN (J * 2)
               K = NELEMS (J)
               IF (NXK (3, K) .EQ. 0) THEN
                  INODE = 2
               ELSEIF (NXK (4, K) .EQ. 0) THEN
                  INODE = 3
               ELSEIF (EIGHT .OR. NINE) THEN
                  INODE = 8
               ELSE
                  INODE = 4
               ENDIF
               IF ( ( (INODE .EQ. 4) .AND.
     &            ( ( (JJ1 .EQ. NXK (1, K)) .AND.
     &            (JJ2 .EQ. NXK (2, K)) ) .OR.
     &            ( (JJ2 .EQ. NXK (1, K)) .AND.
     &            (JJ1 .EQ. NXK (2, K)) ) ) ) .OR.
     &            ( (INODE .EQ. 8) .AND.
     &            ( ( (JJ1 .EQ. NXK (1, K)) .AND.
     &            (JJ2 .EQ. NXK (2, K)) ) .OR.
     &            ( (JJ2 .EQ. NXK (1, K)) .AND.
     &            (JJ1 .EQ. NXK (2, K)) ) ) ) ) THEN
C     &            ( (JJ1 .EQ. NXK (2, K)) .AND.
C     &            (JJ2 .EQ. NXK (3, K)) ) .OR.
C     &            ( (JJ2 .EQ. NXK (2, K)) .AND.
C     &            (JJ1 .EQ. NXK (3, K)) ) ) ) ) THEN
C
                  IF (.NOT. FOUND) THEN
                     FOUND = .TRUE.
                  ENDIF
C
                  JHOLD = JHOLD + 1
                  IHOLD (JHOLD) = MAPGXD (K)
                  IF (JHOLD .EQ. 9) THEN
                     DO 160 II = 1, JHOLD
                        WRITE (IUNIT, 10110, ERR = 290) IHOLD(II), 1,
     &                     PMAG
  160                CONTINUE
                     JHOLD = 0
                  ENDIF
               ENDIF
  170       CONTINUE
            IF (JHOLD .GT. 0) THEN
               DO 180 II = 1, JHOLD
                  WRITE (IUNIT, 10110, ERR = 290) IHOLD(II), 1, PMAG
  180          CONTINUE
            ENDIF
C
C  WRITE OUT THE SIDE 2 ELEMENTS
C
            FOUND = .FALSE.
            JHOLD = 0
            DO 200 J = J1, J2
               JJ1 = NSIDEN ( (J * 2) - 1)
               JJ2 = NSIDEN (J * 2)
               K = NELEMS (J)
               IF (NXK (3, K) .EQ. 0) THEN
                  INODE = 2
               ELSEIF (NXK (4, K) .EQ. 0) THEN
                  INODE = 3
               ELSEIF (EIGHT .OR. NINE) THEN
                  INODE = 8
               ELSE
                  INODE = 4
               ENDIF
               IF ( ( (INODE .EQ. 4) .AND.
     &            ( ( (JJ1 .EQ. NXK (2, K)) .AND.
     &            (JJ2 .EQ. NXK (3, K)) ) .OR.
     &            ( (JJ2 .EQ. NXK (2, K)) .AND.
     &            (JJ1 .EQ. NXK (3, K)) ) ) ) .OR.
     &            ( (INODE .EQ. 8) .AND.
     &            ( ( (JJ1 .EQ. NXK (3, K)) .AND.
     &            (JJ2 .EQ. NXK (4, K)) ) .OR.
     &            ( (JJ2 .EQ. NXK (3, K)) .AND.
     &            (JJ1 .EQ. NXK (4, K)) ) ) ) ) THEN
C     &            ( (JJ1 .EQ. NXK (4, K)) .AND.
C     &            (JJ2 .EQ. NXK (5, K)) ) .OR.
C     &            ( (JJ2 .EQ. NXK (4, K)) .AND.
C     &            (JJ1 .EQ. NXK (5, K)) ) ) ) ) THEN
C
                  IF (.NOT. FOUND) THEN
                     FOUND = .TRUE.
                  ENDIF
C
                  JHOLD = JHOLD + 1
                  IHOLD (JHOLD) = MAPGXD (K)
                  IF (JHOLD .EQ. 9) THEN
                     DO 190 II = 2, JHOLD
                        WRITE (IUNIT, 10110, ERR = 290) IHOLD(II), 2,
     &                     PMAG
  190                CONTINUE
                     JHOLD = 0
                  ENDIF
               ENDIF
  200       CONTINUE
            IF (JHOLD .GT. 0) THEN
               DO 210 II = 1, JHOLD
                  WRITE (IUNIT, 10110, ERR = 290) IHOLD(II), 2, PMAG
  210          CONTINUE
            ENDIF
C
C  WRITE OUT THE SIDE 3 ELEMENTS
C
            FOUND = .FALSE.
            JHOLD = 0
            DO 230 J = J1, J2
               JJ1 = NSIDEN ( (J * 2) - 1)
               JJ2 = NSIDEN (J * 2)
               K = NELEMS (J)
               IF (NXK (3, K) .EQ. 0) THEN
                  INODE = 2
               ELSEIF (NXK (4, K) .EQ. 0) THEN
                  INODE = 3
               ELSEIF (EIGHT .OR. NINE) THEN
                  INODE = 8
               ELSE
                  INODE = 4
               ENDIF
               IF ( ( (INODE .EQ. 4) .AND.
     &            ( ( (JJ1 .EQ. NXK (3, K)) .AND.
     &            (JJ2 .EQ. NXK (4, K)) ) .OR.
     &            ( (JJ2 .EQ. NXK (3, K)) .AND.
     &            (JJ1 .EQ. NXK (4, K)) ) ) ) .OR.
     &            ( (INODE .EQ. 8) .AND.
     &            ( ( (JJ1 .EQ. NXK (5, K)) .AND.
     &            (JJ2 .EQ. NXK (6, K)) ) .OR.
     &            ( (JJ2 .EQ. NXK (5, K)) .AND.
     &            (JJ1 .EQ. NXK (6, K)) ) ) ) ) THEN
C     &            ( (JJ1 .EQ. NXK (6, K)) .AND.
C     &            (JJ2 .EQ. NXK (7, K)) ) .OR.
C     &            ( (JJ2 .EQ. NXK (6, K)) .AND.
C     &            (JJ1 .EQ. NXK (7, K)) ) ) ) ) THEN
C
                  IF (.NOT. FOUND) THEN
                     FOUND = .TRUE.
                  ENDIF
C
                  JHOLD = JHOLD + 1
                  IHOLD (JHOLD) = MAPGXD (K)
                  IF (JHOLD .EQ. 9) THEN
                     DO 220 II = 1, JHOLD
                        WRITE (IUNIT, 10110, ERR = 290) IHOLD(II), 3,
     &                     PMAG
  220                CONTINUE
                     JHOLD = 0
                  ENDIF
               ENDIF
  230       CONTINUE
            IF (JHOLD .GT. 0) THEN
               DO 240 II = 1, JHOLD
                  WRITE (IUNIT, 10110, ERR = 290) IHOLD(II), 3, PMAG
  240          CONTINUE
            ENDIF
C
C  WRITE OUT THE SIDE 4 ELEMENTS
C
            FOUND = .FALSE.
            JHOLD = 0
            DO 260 J = J1, J2
               JJ1 = NSIDEN ( (J * 2) - 1)
               JJ2 = NSIDEN (J * 2)
               K = NELEMS (J)
               IF (NXK (3, K) .EQ. 0) THEN
                  INODE = 2
               ELSEIF (NXK (4, K) .EQ. 0) THEN
                  INODE = 3
               ELSEIF (EIGHT .OR. NINE) THEN
                  INODE = 8
               ELSE
                  INODE = 4
               ENDIF
               IF ( ( (INODE .EQ. 4) .AND.
     &            ( ( (JJ1 .EQ. NXK (4, K)) .AND.
     &            (JJ2 .EQ. NXK (1, K)) ) .OR.
     &            ( (JJ2 .EQ. NXK (4, K)) .AND.
     &            (JJ1 .EQ. NXK (1, K)) ) ) ) .OR.
     &            ( (INODE .EQ. 8) .AND.
     &            ( ( (JJ1 .EQ. NXK (7, K)) .AND.
     &            (JJ2 .EQ. NXK (8, K)) ) .OR.
     &            ( (JJ2 .EQ. NXK (7, K)) .AND.
     &            (JJ1 .EQ. NXK (8, K)) ) ) ) ) THEN
C     &            ( (JJ1 .EQ. NXK (8, K)) .AND.
C     &            (JJ2 .EQ. NXK (1, K)) ) .OR.
C     &            ( (JJ2 .EQ. NXK (8, K)) .AND.
C     &            (JJ1 .EQ. NXK (1, K)) ) ) ) ) THEN
C
                  IF (.NOT. FOUND) THEN
                     FOUND = .TRUE.
                  ENDIF
C
                  JHOLD = JHOLD + 1
                  IHOLD (JHOLD) = MAPGXD (K)
                  IF (JHOLD .EQ. 9) THEN
                     DO 250 II = 1, JHOLD
                        WRITE (IUNIT, 10110, ERR = 290) IHOLD(II), 4,
     &                     PMAG
  250                CONTINUE
                     JHOLD = 0
                  ENDIF
               ENDIF
  260       CONTINUE
            IF (JHOLD .GT. 0) THEN
               DO 270 II = 1, JHOLD
                  WRITE (IUNIT, 10110, ERR = 290) IHOLD(II), 4, PMAG
  270          CONTINUE
            ENDIF
C
  280    CONTINUE
      ENDIF
      WRITE (IUNIT, 10140, ERR = 290)
      CALL MESAGE ('JOE''S ERROR OUTPUT FILE SUCCESSFULLY WRITTEN')
      ERR = .FALSE.
      RETURN
C
C  ERR DURING WRITE PROBLEMS
C
  290 CONTINUE
      CALL MESAGE ('ERR DURING WRITE TO ABAQUS OUTPUT FILE')
      CALL MESAGE ('         - NO FILE SAVED -            ')
      RETURN
C
C
10000 FORMAT (A72)
10010 FORMAT ('*PROSTRFAC 0.0',/,
     &   '*PSTRESS',/,
     &   '*MATERIAL',/,
     &   '*ELASTIC',/,
     &   '30.E+6 .3',/,
     &   '*MATERIAL',/,
     &   '*ELASTIC',/,
     &   '30.E+6 .3')
10020 FORMAT ('*NODE')
10030 FORMAT (I10, ',', 2 (E14.7, ','), E14.7)
10040 FORMAT ('*ELEMENT')
10050 FORMAT (8 (I8, ','), I8)
10060 FORMAT (4 (I10, ','), I10)
10070 FORMAT ('*FIXED')
10080 FORMAT (/,' FOR NODE BOUNDARY FLAG',I5)
10090 FORMAT (/,' FOR ELEMENT BOUNDARY FLAG',I5)
10100 FORMAT (I8, ',', I8)
10110 FORMAT ('EL,', I8, ', P', I1, ',', F10.6)
10120 FORMAT ('*STEP', /, '1.,1.', /, '*POST', /,
     &   '*MAXIT  2', /, 'PRINT  50')
10130 FORMAT ('*PTOL  1.E-2')
10140 FORMAT ('*END')
C
      END
