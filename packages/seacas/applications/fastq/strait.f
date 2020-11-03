C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE STRAIT (MP, ML, MCOM, ICOM, JCOM, CIN, RIN, IIN, KIN,
     &   IDUMP, N, COOR, LCON, LINKP, LINKL)
C***********************************************************************

C  SUBROUTINE STRAIT = STRAIGHTENS LINES IN THE X OR Y DIRECTION

C***********************************************************************

C  SUBROUTINE CALLED BY:
C     FASTQ = A PROGRAM TO QUICKLY PREPARE QMESH INPUT

C***********************************************************************

C  SUBROUTINES CALLED:
C     CHECK  = CHECKS 2 VALUES FOR BEING OUT OF PRESCRIBED BOUNDS

C***********************************************************************

C  VARIABLES USED:
C     IANS   = LOGICAL RESPONSE FROM YES-NO QUESTION
C     ANS    = CHARACTER RESPONSE FOR MENU CHOICE

C***********************************************************************

      DIMENSION COOR (2, MP), LCON (3, ML), LINKP (2, MP)
      DIMENSION LINKL (2, ML), N (29)
      DIMENSION KIN (MCOM), IIN (MCOM), RIN (MCOM)

      CHARACTER*72 CIN (MCOM)
      LOGICAL ADDLNK

      IZ=0
      ADDLNK=.FALSE.

  100 CONTINUE
      IF (N (2).GT.0)THEN
         CALL MESAGE (' ')
         CALL MESAGE ('STRAIGHTEN LINES <I1> THROUGH <I2> PARALLEL TO'//
     &      ' THE <X OR Y> AXIS')
         IF (ICOM.GT.JCOM)THEN
            CALL FREFLD (IZ, IZ, '>', MCOM, IOSTAT, JCOM, KIN, CIN, IIN,
     &         RIN)
            ICOM=1
         ENDIF
         CALL GETI12 (MCOM, ICOM, JCOM, CIN, IIN, KIN, I1, I2, IFOUND)
         IF (IFOUND.GT.0)THEN
            CALL CHECK (I1, I2, N (19))
         ELSE
            RETURN
         ENDIF

C  STRAIGHTEN THE LINE IN THE Y DIRECTION

         IF ( (CIN (ICOM) (1:1) .EQ. 'Y') .OR.
     &      (CIN (ICOM) (1:1) .EQ. 'y'))THEN
            ICOM=ICOM+1
            DO 120 I=I1, I2
               CALL LTSORT (ML, LINKL, I, J, ADDLNK)
               IF (J.GT.0)THEN
                  CALL LTSORT (MP, LINKP, LCON (1, J), IP1, ADDLNK)
                  CALL LTSORT (MP, LINKP, LCON (2, J), IP2, ADDLNK)
                  IF ( (IP1.GT.0) .AND. (IP2.GT.0) .AND. (IP1.LE.N (18))
     &               .AND. (IP2.LE.N (18)) ) THEN
                     WRITE (*, 10000)I, LCON (1, J), COOR (1, IP1),
     &                  COOR (2, IP1),  LCON (2, J), COOR (1, IP2),
     &                  COOR (2, IP2)
  110                CONTINUE
                     IF (ICOM.GT.JCOM)THEN
                        CALL FREFLD (IZ, IZ, 'NEW CONSTANT X VALUE: ',
     &                     MCOM, IOSTAT, JCOM, KIN, CIN, IIN, RIN)
                        ICOM=1
                     ENDIF
                     IF (KIN (ICOM).GT.0)THEN
                        COOR (1, IP1)=RIN (ICOM)
                        COOR (1, IP2)=RIN (ICOM)
                        ICOM=ICOM+1
                     ELSE
                        CALL MESAGE ('BAD INPUT - TRY AGAIN')
                        ICOM=ICOM+1
                        GOTO 110
                     ENDIF
                  ENDIF
               ENDIF
  120       CONTINUE

C  STRAIGHTEN THE LINE IN THE X DIRECTION

         ELSEIF ( (CIN (ICOM) (1:1).EQ.'X') .OR.
     &      (CIN (ICOM) (1:1).EQ.'x')) THEN
            ICOM=ICOM+1
            DO 140 I=I1, I2
               CALL LTSORT (ML, LINKL, I, J, ADDLNK)
               IF (J.GT.0)THEN
                  CALL LTSORT (MP, LINKP, LCON (1, J), IP1, ADDLNK)
                  CALL LTSORT (MP, LINKP, LCON (2, J), IP2, ADDLNK)
                  IF ( (IP1.GT.0) .AND. (IP2.GT.0) .AND. (IP1.LE.N (18))
     &               .AND. (IP2.LE.N (18)) ) THEN
                     WRITE (*, 10000)I, LCON (1, J), COOR (1, IP1),
     &                  COOR (2, IP1), LCON (2, J), COOR (1, IP2),
     &                  COOR (2, IP2)
  130                CONTINUE
                     IF (ICOM.GT.JCOM)THEN
                        CALL FREFLD (IZ, IZ, 'NEW CONSTANT Y VALUE: ',
     &                     MCOM, IOSTAT, JCOM, KIN, CIN, IIN, RIN)
                        ICOM=1
                     ENDIF
                     IF (KIN (ICOM).GT.0)THEN
                        COOR (2, IP1)=RIN (ICOM)
                        COOR (2, IP2)=RIN (ICOM)
                        ICOM=ICOM+1
                     ELSE
                        CALL MESAGE ('BAD INPUT - TRY AGAIN')
                        ICOM=ICOM+1
                        GOTO 130
                     ENDIF
                  ENDIF
               ENDIF
  140       CONTINUE
         ELSE
            CALL MESAGE ('DIRECTION MUST BE "X" OR "Y"')
            ICOM=ICOM+1
         ENDIF
      ELSE
         CALL MESAGE (' ')
         CALL MESAGE ('*------------------------------*')
         CALL MESAGE ('NO LINES IN THE CURRENT DATABASE')
         CALL MESAGE ('*------------------------------*')
         RETURN
      ENDIF
      GOTO 100

10000 FORMAT (' LINE', I5, ' HAS THE FOLLOWING END POINTS:', /,
     &   ' POINT:', I5, '   X COORDINATE', G14.7,  '   Y COORDINATE ',
     &   G14.7, /,
     &   ' POINT:', I5, '   X COORDINATE', G14.7,  '   Y COORDINATE ',
     &   G14.7)

      END
