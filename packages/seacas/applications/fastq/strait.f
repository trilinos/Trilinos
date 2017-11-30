C    Copyright(C) 2014-2017 National Technology & Engineering Solutions of
C    Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
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
C    * Neither the name of NTESS nor the names of its
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

C $Id: strait.f,v 1.2 1998/07/14 18:20:03 gdsjaar Exp $
C $Log: strait.f,v $
C Revision 1.2  1998/07/14 18:20:03  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:16:40  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:16:39  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]STRAIT.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE STRAIT (MP, ML, MCOM, ICOM, JCOM, CIN, RIN, IIN, KIN,
     &   IDUMP, N, COOR, LCON, LINKP, LINKL)
C***********************************************************************
C
C  SUBROUTINE STRAIT = STRAIGHTENS LINES IN THE X OR Y DIRECTION
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     FASTQ = A PROGRAM TO QUICKLY PREPARE QMESH INPUT
C
C***********************************************************************
C
C  SUBROUTINES CALLED:
C     CHECK  = CHECKS 2 VALUES FOR BEING OUT OF PRESCRIBED BOUNDS
C
C***********************************************************************
C
C  VARIABLES USED:
C     IANS   = LOGICAL RESPONSE FROM YES-NO QUESTION
C     ANS    = CHARACTER RESPONSE FOR MENU CHOICE
C
C***********************************************************************
C
      DIMENSION COOR (2, MP), LCON (3, ML), LINKP (2, MP)
      DIMENSION LINKL (2, ML), N (29)
      DIMENSION KIN (MCOM), IIN (MCOM), RIN (MCOM)
C
      CHARACTER*72 CIN (MCOM)
      LOGICAL ADDLNK
C
      IZ=0
      ADDLNK=.FALSE.
C
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
C
C  STRAIGHTEN THE LINE IN THE Y DIRECTION
C
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
C
C  STRAIGHTEN THE LINE IN THE X DIRECTION
C
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
C
10000 FORMAT (' LINE', I5, ' HAS THE FOLLOWING END POINTS:', /,
     &   ' POINT:', I5, '   X COORDINATE', G14.7,  '   Y COORDINATE ',
     &   G14.7, /,
     &   ' POINT:', I5, '   X COORDINATE', G14.7,  '   Y COORDINATE ',
     &   G14.7)
C
      END
