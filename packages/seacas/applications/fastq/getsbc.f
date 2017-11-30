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

C $Id: getsbc.f,v 1.3 2007/07/24 13:10:18 gdsjaar Exp $
C $Log: getsbc.f,v $
C Revision 1.3  2007/07/24 13:10:18  gdsjaar
C Fix problem with boundary condition memory overwrite.
C
C Remove old ls5 and r25 terminal tests
C
C Revision 1.2  1998/11/24 20:45:07  gdsjaar
C Added code to avoid array bound read errors and uninitialized
C variables. In some cases, the correct fix was difficult to determine,
C so added something that looked like it made sense...
C
C This fixes problems with very slow run times on g77-compiled code. It
C was taking an uninitialized variable to be INT_MAX instead of zero
C which resulted in lots of iterations through a loop. This variable was
C initialized to zero since that is what it was being set to on the sun
C and when compiled with fort77 (f2c-based).  Gives the exact same mesh
C on linux and sun for several test cases.
C
C Revision 1.1.1.1  1990/11/30 11:08:44  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:08:42  gdsjaar
c Initial revision
c 
C
CC* FILE: [.QMESH]GETSBC.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GETSBC (MXND, MXNPER, NPER, NL, ML, MAXSBC, MAXPRM,
     &   NPRM, NID, LISTL, XN, YN, NUID, LXK, KXL, NXL, LSTSBC, NPERIM,
     &   KSBC, LCON, ISBOUN, LINKL, NSPF, IFSB, LISTSB, LINKSB, LLL,
     &   BAR, ERR)
C***********************************************************************
C
C  SUBROUTINE GETSBC = GETS THE SIDE BOUNDARY LIST
C
C***********************************************************************
C
      DIMENSION NID (MXNPER, MAXPRM), NPERIM (MAXPRM)
      DIMENSION LISTL (NL), XN (MXND), YN (MXND), NUID (MXND)
      DIMENSION LXK (4, MXND), KXL (2, MXND*3), NXL (2, MXND*3)
      DIMENSION LCON (3, ML), ISBOUN (ML), LINKL (2, ML)
      DIMENSION NSPF (ML), IFSB (ML), LISTSB (2, ML), LINKSB (2, ML)
      DIMENSION LSTSBC (MAXSBC), NODES (4)
C
      LOGICAL EXTER, ERR, CCW, BAR, ADDLNK
C
      ERR = .TRUE.
      CCW = .TRUE.
      ADDLNK = .FALSE.
      NPERIM (1) = NPER
C
      DO 110 II = 1, NPRM
         DO 100 I = 1, NPERIM(II)
            IF (BAR) THEN
               NID (I, II) = IABS (NUID (I))
            ELSE
               IF (NID (I, II) .LT. 0) NID (I, II) = - NID (I, II)
            ENDIF
  100    CONTINUE
  110 CONTINUE
C
C  SORT THROUGH AND PICK OFF ELEMENTS WITH SIDE BOUNDARY CONDITIONS
C
      DO 240 I = 1, LLL
         IF (BAR) THEN
            I1 = LXK (1, I)
            I2 = LXK (2, I)
         ELSE
            I1 = NXL (1, I)
            I2 = NXL (2, I)
         ENDIF
C
C  SEE IF THE LINE IS CLEARLY INTERIOR
C
         IF (I1 .GT. 0 .AND. I2 .GT. 0) THEN
           if ((NUID (I1) .NE. 0) .AND. (NUID (I2) .NE. 0)) THEN
            LTEST = 0
            EXTER = .FALSE.
C
C  CHECK AGAINST THE PERIMETER LIST TO SEE IF IT IS TRULLY EXTERIOR
C
            DO 130 JJ  =  1, NPRM
               DO 120 J = 1, NPERIM (JJ)
                  IF (ABS (NUID (I1)) .EQ. NID (J, JJ)) THEN
                     IF (J .EQ. 1) THEN
                        J1 = J + 1
                        J2 = NPERIM(JJ)
                     ELSEIF (J .EQ. NPERIM(JJ)) THEN
                        J1 = J - 1
                        J2 = 1
                     ELSE
                        J1 = J - 1
                        J2 = J + 1
                     ENDIF
                     IF ( (ABS (NUID (I2)) .EQ. NID (J1, JJ)) .OR.
     &                  (ABS (NUID (I2)) .EQ. NID (J2, JJ)) )
     &                  EXTER = .TRUE.
                     GOTO 140
                  ENDIF
  120          CONTINUE
  130       CONTINUE
  140       CONTINUE
            IF (EXTER) THEN
C
C  FIND THE LINE NUMBER IT BELONGS TO
C
               IF (ABS (NUID (I1)) .GT. 1000000000) THEN
                  LTEST =  (ABS (NUID (I1)) - 1000000000) / 100000
               ELSEIF (ABS (NUID (I2)) .GT. 1000000000) THEN
                  LTEST =  (ABS (NUID (I2)) - 1000000000) / 100000
               ELSE
                  NSUM = ABS (NUID (I1)) + ABS (NUID (I2))
                  DO 150 J = 1, NL
                     CALL LTSORT (ML, LINKL, LISTL (J), K, ADDLNK)
                     IF ((LCON (1, K) + LCON (2, K)) .EQ. NSUM) THEN
                        IF (( (LCON (1, K) .EQ. ABS (NUID (I1))) .AND.
     +                     (LCON (2, K) .EQ. ABS (NUID (I2)))) .OR.
     +                     ((LCON (1, K) .EQ. ABS (NUID (I2))) .AND.
     +                     (LCON (2, K) .EQ. ABS (NUID (I1))))) THEN
                           LTEST = LISTL (J)
                           GOTO 160
                        ENDIF
                     ENDIF
  150             CONTINUE
  160             CONTINUE
               ENDIF
C
C  FIND THE ELEMENT BOUNDARY FLAG IF THERE IS ONE
C
               IF (LTEST.LE.0) THEN
                  CALL MESAGE (' ERROR IN SEARCHING NXL FOR '//
     &               'ELEMENT BCC')
                  RETURN
               ELSE
                  CALL LTSORT (ML, LINKL, LTEST, J, ADDLNK)
                  IF (ISBOUN (J) .GT. 0) THEN
                     IFLAG = ISBOUN (J)
C
C  CHECK TO MAKE SURE LINE IS LINKED TO FLAG
C  AND GET THE NEXT LINK  (NFLAG)
C
                     CALL LTSORT (ML, LINKSB, IFLAG, L, ADDLNK)
                     DO 170 JJ = IFSB (L), IFSB (L) + NSPF (L) - 1
                        IF (LISTSB (1, JJ) .LT. 0) THEN
                           CALL MESAGE ('PROBLEMS WITH SIDES IN '//
     &                        'FLAG LIST IN GETSBC')
                        ELSE
                           IF (LISTSB (1, JJ) .EQ. LTEST) THEN
                              NFLAG = LISTSB (2, JJ)
                              GOTO 180
                           ENDIF
                        ENDIF
  170                CONTINUE
                     WRITE (*, 10000)IFLAG
                     RETURN
  180                CONTINUE
                     IF (BAR) THEN
                        NELEM = I
                     ELSE
                        NELEM = KXL (1, I)
                        IF (NELEM .EQ. 0)NELEM = KXL (2, I)
                     ENDIF
                     KSBC = KSBC + 1
                     LSTSBC (KSBC) = - IFLAG
                     KSBC = KSBC + 1
                     if (ksbc .gt. maxsbc) stop 'maxsbc error'
                     LSTSBC (KSBC) = NELEM
C
C  GET THE CORRECT ELEMENT SIDE
C
                     IF (BAR) THEN
                        JSIDE = 1
                     ELSE
                        CALL GNXKA (MXND, XN, YN, NELEM, NODES, AREA,
     &                     LXK, NXL, CCW)
                        DO 190 J = 1, 4
                           IF (I1 .EQ. NODES (J)) THEN
                              JP1 = J + 1
                              JM1 = J - 1
                              IF (JP1 .EQ. 5)JP1 = 1
                              IF (JM1 .EQ. 0)JM1 = 4
                              IF (I2 .EQ. NODES (JP1)) THEN
                                 JSIDE = J
                                 GOTO 200
                              ELSEIF (I2 .EQ. NODES (JM1)) THEN
                                 JSIDE = JM1
                                 GOTO 200
                              ENDIF
                           ENDIF
  190                   CONTINUE
                        WRITE (*, 10010)NELEM
                        RETURN
  200                   CONTINUE
                     ENDIF
                     KSBC = KSBC + 1
                     LSTSBC (KSBC) = JSIDE
C
C  SEE IF ANY MORE FLAGS ARE ATTACHED TO THIS SIDE
C
  210                CONTINUE
                     IF (NFLAG .GT. 0) THEN
C
C  CHECK TO MAKE SURE LINE IS LINKED TO FLAG
C  AND GET THE NEXT LINK  (NFLAG)
C
                        IFLAG = NFLAG
                        CALL LTSORT (ML, LINKSB, IFLAG, L, ADDLNK)
                        DO 220 JJ = IFSB (L), IFSB (L) + NSPF (L)
                           IF (LISTSB (1, JJ) .EQ. LTEST) THEN
                              NFLAG = LISTSB (2, JJ)
                              GOTO 230
                           ENDIF
  220                   CONTINUE
                        WRITE (*, 10000)IFLAG
                        RETURN
  230                   CONTINUE
                        KSBC = KSBC + 1
                        LSTSBC (KSBC) = - IFLAG
                        KSBC = KSBC + 1
                        LSTSBC (KSBC) = NELEM
                        KSBC = KSBC + 1
                        LSTSBC (KSBC) = JSIDE
                        GOTO 210
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
       END IF
  240 CONTINUE
C
      ERR = .FALSE.
      RETURN
C
10000 FORMAT (' SIDE BOUNDARY FLAG', I5, ' IS NOT PROPERLY LINKED')
10010 FORMAT (' ERROR FINDING CORRECT BOUNDARY SIDE ON ELEMENT', I5)
      END
