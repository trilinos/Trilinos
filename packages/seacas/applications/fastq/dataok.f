C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Governement retains certain rights in this software.
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

C $Id: dataok.f,v 1.2 1991/05/10 17:40:36 gdsjaar Exp $
C $Log: dataok.f,v $
C Revision 1.2  1991/05/10 17:40:36  gdsjaar
C Changed VMS JNINT to ANSI NINT, but then had
C to change variable NINT to KNINT
C
c Revision 1.1.1.1  1990/11/30  11:05:44  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:05:42  gdsjaar
c Initial revision
c 
C
CC* FILE: [.QMESH]DATAOK.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE DATAOK (MP, ML, MS, MR, L, KNUM, COOR, ILINE, LTYPE,
     &   KNINT, LCON, NLPS, IFLINE, ILLIST, NSPR, IFSIDE, ISLIST, LINKP,
     &   LINKL, LINKS, SIZE, ERRCHK, ERR)
C***********************************************************************
C
C  SUBROUTINE FILLOK = CHECKS TO MAKE SURE NONEXISTANT DATA IS NOT
C                      BEING REFERENCED IN THE REGION DEFINITIONS
C
C***********************************************************************
C
      DIMENSION COOR (2, MP), LINKP (2, MP)
      DIMENSION ILINE (ML), LTYPE (ML), KNINT (ML), LCON (3, ML)
      DIMENSION LINKL (2, ML)
      DIMENSION NLPS (MS), IFLINE (MS), ILLIST (MS*3), LINKS (2, MS)
      DIMENSION NSPR (MR), IFSIDE (MR), ISLIST (MR*4)
C
      LOGICAL ERR, ADDLNK, ERRCHK
C
      ERR = .TRUE.
      ADDLNK = .FALSE.
C
      DO 130 I = IFSIDE (L), IFSIDE (L) + NSPR (L)-1
C
C  CHECK TO MAKE SURE REGION'S SIDE DEFINITIONS ARE ALL THERE
C
         IF (ISLIST (I).GT.0)THEN
            II = ISLIST (I)
            CALL LTSORT (MS, LINKS, II, IPNTR, ADDLNK)
            IF (IPNTR.LE.0) THEN
               IF (ERRCHK) THEN
                  WRITE (*, 10000)KNUM, II
                  RETURN
               ELSE
                  GOTO 120
               ENDIF
            END IF
C
C  CHECK TO MAKE SURE SIDE'S LINE DEFINITIONS ARE ALL THERE
C
            CALL LTSORT (MS, LINKS, II, JJ, ADDLNK)
            DO 110 J = IFLINE (JJ), IFLINE (JJ) + NLPS (JJ)-1
               KK = ILLIST (J)
               CALL LTSORT (ML, LINKL, KK, LL, ADDLNK)
               IF ((KK.LE.0) .OR. (LL.LE.0)) THEN
                  IF (ERRCHK) THEN
                     WRITE (*, 10010)II, KK
                     RETURN
                  ELSE
                     GOTO 100
                  ENDIF
               END IF
C
C  CHECK TO MAKE SURE LINE'S POINT DEFINITIONS ARE ALL THERE
C
               I1 = LCON (1, LL)
               I2 = LCON (2, LL)
               I3 = LCON (3, LL)
               CALL LTSORT (MP, LINKP, I1, J1, ADDLNK)
               CALL LTSORT (MP, LINKP, I2, J2, ADDLNK)
               IF (I3 .NE. 0) THEN
                  CALL LTSORT (MP, LINKP, IABS (I3), J3, ADDLNK)
               ELSE
                  J3 = 0
               END IF
C
               IF ((I1.LE.0) .OR. (J1.LE.0)) THEN
                  IF (ERRCHK) THEN
                     WRITE (*, 10030)KK, I1
                     RETURN
                  ELSE
                     GOTO 100
                  ENDIF
               ELSEIF ((I2.LE.0) .OR. (J2.LE.0)) THEN
                  IF (ERRCHK) THEN
                     WRITE (*, 10030)KK, I2
                     RETURN
                  ELSE
                     GOTO 100
                  ENDIF
               ELSEIF ((LTYPE (LL) .NE. 1) .AND. ((I3 .EQ. 0) .OR.
     &            (J3.LE.0))) THEN
                  IF (ERRCHK) THEN
                     WRITE (*, 10030)KK, I3
                     RETURN
                  ELSE
                     GOTO 100
                  ENDIF
               END IF
C
C  CHECK TO INSURE AN INTERAL ASSIGNMENT
C
               IF (IABS (KNINT (LL)) .EQ. 0) THEN
                  IF (I3 .LT. 0)J3 = -J3
                  CALL LINLEN (MP, COOR, LINKP, KNUM, ILINE(LL),
     &               LTYPE(LL), I3, J1, J2, J3, DIST, ERR)
                  IF (ERR) THEN
                     IF (ERRCHK) THEN
                        WRITE (*, 10020)KK, IABS (KNINT (LL))
                        RETURN
                     ELSE
                        GOTO 100
                     ENDIF
                  ELSE
                     IF (SIZE .LE. 0.) THEN
                        KNINT (LL) = 1
                     ELSE
                        KNINT (LL) = MAX0 (NINT (DIST/SIZE), 1)
                     END IF
                  END IF
               END IF
  100          CONTINUE
  110       CONTINUE
C
C  CHECK TO MAKE SURE REGION'S LINE DEFINITIONS ARE ALL THERE
C
         ELSEIF (ISLIST (I) .LT. 0) THEN
            KK = IABS (ISLIST (I))
            CALL LTSORT (ML, LINKL, KK, LL, ADDLNK)
            IF ( (KK .LE. 0)  .OR. (LL .LE. 0) ) THEN
               IF (ERRCHK) THEN
                  WRITE (*, 10010)KNUM, KK
                  RETURN
               ELSE
                  GOTO 120
               ENDIF
            END IF
C
C  CHECK TO MAKE SURE LINE'S POINT DEFINITIONS ARE ALL THERE
C
            I1 = LCON (1, LL)
            I2 = LCON (2, LL)
            I3 = LCON (3, LL)
            CALL LTSORT (MP, LINKP, I1, J1, ADDLNK)
            CALL LTSORT (MP, LINKP, I2, J2, ADDLNK)
            IF (I3 .NE. 0) THEN
               CALL LTSORT (MP, LINKP, IABS (I3), J3, ADDLNK)
            ELSE
               J3 = 0
            END IF
C
            IF ((I1.LE.0) .OR. (J1.LE.0)) THEN
               IF (ERRCHK) THEN
                  WRITE (*, 10030)KK, I1
                  RETURN
               ELSE
                  GOTO 120
               ENDIF
            ELSEIF ((I2.LE.0) .OR. (J2.LE.0)) THEN
               IF (ERRCHK) THEN
                  WRITE (*, 10030)KK, I2
                  RETURN
               ELSE
                  GOTO 120
               ENDIF
            ELSEIF ((LTYPE (LL) .NE. 1) .AND. ((I3 .EQ. 0) .OR.
     &         (J3.LE.0))) THEN
               IF (ERRCHK) THEN
                  WRITE (*, 10030)KK, I3
                  RETURN
               ELSE
                  GOTO 120
               ENDIF
            END IF
C
C  CHECK TO MAKE SURE INTERVAL ASSIGNMENT IS HANDLED
C
            IF (IABS (KNINT (LL)) .EQ. 0) THEN
C
C**MBS/29-JUN-1989/ DO NOT NEGATE POINTER TO CENTER OF CLOCKWISE ARC
C              IF (I3 .LT. 0)J3 = -J3
               CALL LINLEN (MP, COOR, LINKP, KNUM, ILINE(LL),
     &            LTYPE(LL), I3, J1, J2, J3, DIST, ERR)
               IF (ERR) THEN
                  IF (ERRCHK) THEN
                     WRITE (*, 10020)KK, IABS (KNINT (LL))
                     RETURN
                  ELSE
                     GOTO 120
                  ENDIF
               ELSE
                  IF (SIZE .LE. 0.) THEN
                     KNINT (LL) = 1
                  ELSE
                     KNINT (LL) = MAX0 (NINT (DIST/SIZE), 1)
                  END IF
               END IF
            END IF
C
C  A ZERO SIDE NUMBER HAS BEEN FOUND
C
         ELSE
            IF (ERRCHK) THEN
               WRITE (*, 10000)KNUM, ISLIST (I)
            ELSE
               GOTO 120
            ENDIF
         END IF
  120    CONTINUE
  130 CONTINUE
C
C  ALL DEFINITIONS ARE IN ORDER
C
      ERR = .FALSE.
      RETURN
C
10000 FORMAT (' FOR REGION:', I5, ' SIDE:', I5, ' DOES NOT EXIST')
10010 FORMAT (' FOR SIDE:', I5, ' LINE:', I5, ' DOES NOT EXIST')
10020 FORMAT (' FOR LINE:', I5, ' INTERVAL OF:', I5, ' IS NOT WORKING')
10030 FORMAT (' FOR LINE:', I5, ' POINT:', I5, ' DOES NOT EXIST')
C
      END
