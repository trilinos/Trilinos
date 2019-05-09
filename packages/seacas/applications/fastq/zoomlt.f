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

C $Id: zoomlt.f,v 1.3 2007/07/24 13:10:18 gdsjaar Exp $
C $Log: zoomlt.f,v $
C Revision 1.3  2007/07/24 13:10:18  gdsjaar
C Fix problem with boundary condition memory overwrite.
C
C Remove old ls5 and r25 terminal tests
C
C Revision 1.2  1998/07/14 18:20:20  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:18:05  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:18:04  gdsjaar
c Initial revision
c
C
CC* FILE: [.MAIN]ZOOMLT.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE ZOOMLT (MCOM, ICOM, JCOM, CIN, RIN, IIN, KIN, IDUMP,
     &   DRAWN, ALPHA, DEV1, X1, X2, Y1, Y2, XX1, XX2, YY1, YY2, XMIN1,
     &   XMAX1, YMIN1, YMAX1, XMIN, XMAX, YMIN, YMAX)
C***********************************************************************
C
C  ZOOMPL = SUBROUTINE TO INPUT NEW ZOOM LIMITS
C
C***********************************************************************
C
      DIMENSION KIN(MCOM), IIN(MCOM), RIN(MCOM)
C
      CHARACTER*72 CIN(MCOM)
      CHARACTER*3 DEV1, ANS
C
      LOGICAL DRAWN, ALPHA
C
      IF ((ICOM .LE. JCOM) .AND. (DRAWN) .AND.
     &   ((CIN(ICOM)(1:1) .EQ. 'C') .OR.
     &   (CIN(ICOM)(1:1) .EQ. 'c')) .AND.
     &   (.NOT.ALPHA)) THEN
         CIN(ICOM) = 'PLOT'
C
C  USE CURSOR INPUT FROM THE SCREEN
C
         CALL MESAGE (' ')
         CALL MESAGE (' ')
         CALL MESAGE ('LOCATE ONE CORNER WITH CURSOR')
         CALL MESAGE ('THEN HIT ANY KEY')
C         X1 = .45
C         Y1 = .325
         CALL PLTCRS (XTEST1, YTEST1, ANS)
C         XDRAW = ABS ( (XTEST1 * (XX2 - XX1))) + XX1
C         YDRAW = ABS ( ((YTEST1 / .75) * (YY2 - YY1)) ) + YY1
         CALL PLTMOV (0., YTEST1)
         CALL PLTDRW (1., YTEST1)
         CALL PLTMOV (XTEST1, 0.)
         CALL PLTDRW (XTEST1, .75)
         CALL PLTFLU
C         X2 = MAX( X1+.05, .55)
C         Y2 = MAX( Y1+.05, .425)
         CALL PLTCRS (XTEST2, YTEST2, ANS)
         CALL PLTMOV (0., YTEST2)
         CALL PLTDRW (1., YTEST2)
         CALL PLTMOV (XTEST2, 0.)
         CALL PLTDRW (XTEST2, .75)
         CALL MESAGE ('LOCATE THE OTHER CORNER WITH CURSOR')
         CALL MESAGE ('THEN HIT ANY KEY')
         CALL PLTFLU
         X1 = MIN (XTEST1, XTEST2)
         X2 = MAX (XTEST1, XTEST2)
         Y1 = MIN (YTEST1, YTEST2)
         Y2 = MAX (YTEST1, YTEST2)
         XMIN = ABS ( (X1 * (XX2 - XX1))) + XX1
         XMAX = ABS ( (X2 * (XX2 - XX1))) + XX1
         YMIN = ABS ( ((Y1 / .75) * (YY2 - YY1)) ) + YY1
         YMAX = ABS ( ((Y2 / .75) * (YY2 - YY1)) ) + YY1
C
C  USE USER INPUT FROM THE KEYPAD
C
      ELSE
         IF ((CIN(ICOM)(1:1) .EQ. 'C') .OR.
     &      (CIN(ICOM)(1:1) .EQ. 'c')) THEN
            ICOM = ICOM+1
            CALL MESAGE (' ')
            CALL MESAGE ('NO CURRENT PLOT FOR CURSOR ZOOM')
            CALL MESAGE ('CURRENT PLOT LIMITS UNCHANGED')
            CALL MESAGE ('* IN OTHER WORDS ... PLOT FIRST (P) '//
     &         'AND THEN ZOOM (Z,C) *')
C
C  SEE IF ANY OF THE VALUES ARE REDEFINED
C
         ELSE IF ( (ICOM .LE. JCOM)  .AND.
     &      ( (KIN(ICOM) .GT. 0)    .OR.  (KIN(ICOM+1) .GT. 0)  .OR.
     &      (KIN(ICOM+2) .GT. 0)  .OR.  (KIN(ICOM+3) .GT. 0) ) ) THEN
            IF (KIN(ICOM) .GT. 0) XMIN = RIN(ICOM)
            ICOM = ICOM+1
            IF (ICOM .LE. JCOM) THEN
               IF (KIN(ICOM) .GT. 0) XMAX = RIN(ICOM)
               ICOM = ICOM+1
               IF (ICOM .LE. JCOM) THEN
                  IF (KIN(ICOM) .GT. 0) YMIN = RIN(ICOM)
                  ICOM = ICOM+1
                  IF (ICOM .LE. JCOM) THEN
                     IF (KIN(ICOM) .GT. 0) YMAX = RIN(ICOM)
                     ICOM = ICOM+1
                  END IF
               END IF
            END IF
         ELSE
            XMIN = XMIN1
            YMIN = YMIN1
            XMAX = XMAX1
            YMAX = YMAX1
            CALL MESAGE (' ')
            CALL MESAGE ('ZOOM LIMITS RESET TO PLOT EXTREMES')
         END IF
      END IF
C
      RETURN
C
      END
