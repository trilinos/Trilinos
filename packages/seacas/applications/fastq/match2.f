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

C $Id: match2.f,v 1.1 1990/11/30 11:11:53 gdsjaar Exp $
C $Log: match2.f,v $
C Revision 1.1  1990/11/30 11:11:53  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]MATCH2.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE MATCH2 (MXND, MLN, XN, YN, NXL, LXN, LNODES, ANGLE,
     &   N0, N1, N2, N3, N0TEST, N1TEST, N2TEST, N3TEST, I1, I2,
     &   J1, J2, KOUNTL, LMATCH, KOUNT2, NODE, U, W, NLOOP, PMATCH, ERR)
C***********************************************************************
C
C  SUBROUTINE MATCH2 = MATCHES UP THE BEST PAIR OF LINES FOR COLLAPSING
C
C***********************************************************************
C
      DIMENSION XN (MXND), YN (MXND), NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION ANGLE (MXND), LNODES (MLN, MXND)
C
      LOGICAL LMATCH, CORNP, SIDEP, BWINS, FWINS, MATCHK, ERR, PMATCH
C
      ERR = .FALSE.
      LMATCH = .TRUE.
      BWINS = .FALSE.
      FWINS = .FALSE.
C
C  MAKE SURE THAT AN ODD NUMBER OR SMALL NUMBER IN A LOOP HAS NOT BEEN
C  CUT OFF, AND IF IT HAS ADJUST THE INTERSECTION ACCORDINGLY.
C
C  FIRST CHECK A 2 NODE LOOP - A HIGHLY UNLIKELY CONDITION
C
      IF (PMATCH) THEN
         I1 = N1
         I2 = N2
         J1 = N1TEST
         J2 = N2TEST
         KOUNTL = KOUNT2 - 1
      ELSEIF (KOUNT2 .EQ. 2) THEN
         IF ((CORNP (ANGLE (N2)) ) .AND. (LXN (4, N2) .NE. 0) ) THEN
            I1 = N1
            I2 = N2
            J1 = N2TEST
            J2 = N3TEST
            KOUNTL = 2
         ELSEIF ((CORNP (ANGLE (N1TEST)) ) .AND.
     &      (LXN (4, N1TEST) .NE. 0) ) THEN
            I1 = N1
            I2 = N2
            J1 = N2TEST
            J2 = N3TEST
            KOUNTL = 2
         ELSE
            IF (CORNP (ANGLE (N2)) ) THEN
               NODE = N2
            ELSEIF (CORNP (ANGLE (N1TEST)) ) THEN
               NODE = N1TEST
            ELSE
               NODE = N2
            ENDIF
            LMATCH = .FALSE.
            GOTO 100
         ENDIF
C
C  NEXT CHECK A 3 NODE LOOP - THIS IS A MUCH MORE PLAUSIBLE CONDITION
C
      ELSEIF (KOUNT2 .EQ. 3) THEN
C
C  CHECK FOR A 3-1 SEMICIRCLE BEING FORMED EITHER WAY
C
         IF ( ( CORNP (ANGLE (N2)) ) .AND.
     &      ( LXN (4, N2) .NE. 0) .AND.
     &      ( SIDEP (ANGLE (N3)) ) ) THEN
            I1 = N0
            I2 = N1
            J1 = N3
            J2 = N1TEST
            KOUNTL = 2
         ELSEIF ( ( CORNP (ANGLE (N1TEST)) ) .AND.
     &      ( LXN (4, N1TEST) .NE. 0) .AND.
     &      ( SIDEP (ANGLE (N3)) ) ) THEN
            I1 = N2
            I2 = N3
            J1 = N2TEST
            J2 = N3TEST
            KOUNTL = 2
C
C  JUST PUT IT AT TWO NODES LEFT
C
         ELSE
            I1 = N1
            I2 = N2
            J1 = N1TEST
            J2 = N2TEST
            KOUNTL = 2
         ENDIF
C
C  NODE LOOP FOR AN EVEN NUMBER OF SPLITS - THE MATCH IS
C  NOT FINE, SO A SHIFT ONE WAY OR THE OTHER IS NEEDED
C
      ELSEIF (MOD (KOUNT2, 2) .EQ. 0) THEN
         I1 = N1
         I2 = N2
C
         XI = XN (I2) - XN (I1)
         YI = YN (I2) - YN (I1)
         XJF = XN (N2TEST) - XN (N3TEST)
         YJF = YN (N2TEST) - YN (N3TEST)
         XJB = XN (N0TEST) - XN (N1TEST)
         YJB = YN (N0TEST) - YN (N1TEST)
C
         FDOT = ( (XI * XJF) + (YI * YJF) ) /
     &      ( SQRT ( (XI * XI) + (YI * YI) ) *
     &      SQRT ( (XJF * XJF) + (YJF * YJF) ) )
         D1F = SQRT ( (XN (N2TEST) - XN (I2)) ** 2 +
     &      (YN (N2TEST) - YN (I2)) ** 2 )
         D2F = SQRT ( (XN (N3TEST) - XN (I1)) ** 2 +
     &      (YN (N3TEST) - YN (I1)) ** 2 )
         DF = (D1F + D2F) * .5
C
         BDOT = ( (XI * XJB) + (YI * YJB) ) /
     &      ( SQRT ( (XI * XI) + (YI * YI) ) *
     &      SQRT ( (XJB * XJB) + (YJB * YJB) ) )
         D1B = SQRT ( (XN (N0TEST) - XN (I2)) ** 2 +
     &      (YN (N0TEST) - YN (I2)) ** 2 )
         D2B = SQRT ( (XN (N1TEST) - XN (I1)) ** 2 +
     &      (YN (N1TEST) - YN (I1)) ** 2 )
         DB = (D1B + D2B) * .5
C
C  NOW COMPARE A FORWARD OR BACKWARD SHIFT AND PICK THE MOST
C  APPROPRIATE ONE BASED ON ANGLE COSINE AND END DISTANCES
C  IF ANY STICK OUT AS THE MOST APPROPRIATE
C
         IF ((FDOT .GT. BDOT) .AND. (DF .LE. DB)) THEN
            J1 = N2TEST
            J2 = N3TEST
            KOUNTL = KOUNT2
         ELSEIF ((BDOT .GT. FDOT) .AND. (DB .LE. DF) .AND.
     &      (KOUNT2 .GT. 4)) THEN
            J1 = N0TEST
            J2 = N1TEST
            KOUNTL = KOUNT2 - 2
         ELSEIF (ABS (ABS( ACOS (BDOT)) - ABS (ACOS (FDOT))) .LE.
     &      .3490659) THEN
            IF ((DF .LE. DB) .OR. (KOUNT2 .LE. 4)) THEN
               J1 = N2TEST
               J2 = N3TEST
               KOUNTL = KOUNT2
            ELSE
               J1 = N0TEST
               J2 = N1TEST
               KOUNTL = KOUNT2 - 2
            ENDIF
C
C  NONE STICK OUT AS THE OVIOUS WINNER - TAKE ONE BASED ON
C  INTERSECTION PORTIONS
C
         ELSE
            IF (U .LT. .5) THEN
               IF ((W .LT. .5) .AND. (KOUNT2 .GT. 4)) THEN
                  J1 = N0TEST
                  J2 = N1TEST
                  KOUNTL = KOUNT2 - 2
               ELSE
                  J1 = N2TEST
                  J2 = N3TEST
                  KOUNTL = KOUNT2
               ENDIF
            ELSE
               IF ((W .LT. .5) .AND. (KOUNT2 .GT. 4)) THEN
                  J1 = N0TEST
                  J2 = N1TEST
                  KOUNTL = KOUNT2 - 2
               ELSE
                  J1 = N2TEST
                  J2 = N3TEST
                  KOUNTL = KOUNT2
               ENDIF
            ENDIF
         ENDIF
C
C  NODE LOOP FOR AN ODD NUMBER OF SPLITS - THE MATCH IS FINE
C
      ELSE
         I1 = N1
         I2 = N2
         J1 = N1TEST
         J2 = N2TEST
         KOUNTL = KOUNT2 - 1
      ENDIF
C
C  NOW THAT THE INITIAL MATCH IS MADE, CHECK MOVING BOTH SIDES
C  FORWARD OR BACKWARD ONE NOTCH AND SEE IF THAT MATCH MAKES MORE SENSE
C  THEN THE CURRENT MATCH
C
      IFOR1 = I2
      IFOR2 = LNODES (3, I2)
      IBAC1 = LNODES (2, I1)
      IBAC2 = I1
C
      JFOR1 = J1
      JFOR2 = LNODES (2, J1)
      JBAC1 = LNODES (3, J2)
      JBAC2 = J2
C
C  NOW CALCULATE THE CROSS PRODUCT AND END DISTANCES
C
      XIF = XN (IFOR2) - XN (IFOR1)
      YIF = YN (IFOR2) - YN (IFOR1)
      XJF = XN (JFOR2) - XN (JFOR1)
      YJF = YN (JFOR2) - YN (JFOR1)
      FDOT = ( (XIF * XJF) + (YIF * YJF) ) /
     &   ( SQRT ( (XIF * XIF) + (YIF * YIF) ) *
     &   SQRT ( (XJF * XJF) + (YJF * YJF) ) )
      D1F = SQRT ( (XN (IFOR1) - XN (JFOR1)) ** 2 +
     &   (YN (IFOR1) - YN (JFOR1)) ** 2 )
      D2F = SQRT ( (XN (IFOR2) - XN (JFOR2)) ** 2 +
     &   (YN (IFOR2) - YN (JFOR2)) ** 2 )
      DF = (D1F + D2F) * .5
C
      XIB = XN (IBAC2) - XN (IBAC1)
      YIB = YN (IBAC2) - YN (IBAC1)
      XJB = XN (JBAC2) - XN (JBAC1)
      YJB = YN (JBAC2) - YN (JBAC1)
      BDOT = ( (XIB * XJB) + (YIB * YJB) ) /
     &   ( SQRT ( (XIB * XIB) + (YIB * YIB) ) *
     &   SQRT ( (XJB * XJB) + (YJB * YJB) ) )
      D1B = SQRT ( (XN (IBAC1) - XN (JBAC1)) ** 2 +
     &   (YN (IBAC1) - YN (JBAC1)) ** 2 )
      D2B = SQRT ( (XN (IBAC2) - XN (JBAC2)) ** 2 +
     &   (YN (IBAC2) - YN (JBAC2)) ** 2 )
      DB = (D1B + D2B) * .5
C
      XI = XN (I2) - XN (I1)
      YI = YN (I2) - YN (I1)
      XJ = XN (J1) - XN (J2)
      YJ = YN (J1) - YN (J2)
      DOT = ( (XI * XJ) + (YI * YJ) ) /
     &   ( SQRT ( (XI * XI) + (YI * YI) ) *
     &   SQRT ( (XJ * XJ) + (YJ * YJ) ) )
      D1 = SQRT ( (XN (I1) - XN (J2)) ** 2 +
     &   (YN (I1) - YN (J2)) ** 2 )
      D2 = SQRT ( (XN (I2) - XN (J1)) ** 2 +
     &   (YN (I2) - YN (J1)) ** 2 )
      D0 = (D1 + D2) * .5
C
C  NOW COMPARE TO SEE IF ANOTHER COMBINATION MAKES BETTER SENSE
C
      IF ( ( ((FDOT .GT. DOT) .AND. (DF .LT. D0)) .OR.
     &   ((.6 * FDOT .GT. DOT) .AND. (DF * .7 .LT. D0)) .OR.
     &   ((.2 * FDOT .GT. DOT) .AND. (DF * .5 .LT. D0)) ) .AND.
     &   (KOUNTL .GT. 4) ) THEN
         FWINS = .TRUE.
         D0 = DF
         DOT = FDOT
      ENDIF
      IF ( ((BDOT .GT. DOT) .AND. (DB .LT. D0)) .OR.
     &   ((.6 * BDOT .GT. DOT) .AND. (DB * .7 .LT. D0)) .OR.
     &   ((.2 * BDOT .GT. DOT) .AND. (DB * .5 .LT. D0)) .AND.
     &   (NLOOP - KOUNTL - 2 .GT. 4) ) THEN
         BWINS = .TRUE.
      ENDIF
C
      IF (BWINS) THEN
         I1 = IBAC1
         I2 = IBAC2
         J1 = JBAC2
         J2 = JBAC1
         KOUNTL = KOUNTL + 2
      ELSEIF (FWINS) THEN
         I1 = IFOR1
         I2 = IFOR2
         J1 = JFOR2
         J2 = JFOR1
         KOUNTL = KOUNTL - 2
      ENDIF
C
C  NOW CHECK THAT TWO BOUNDARY LINES OR LINES CONNECTED TO THE
C  BOUNDARY ARE NOT BEING JOINED INAPPROPRIATELY
C
      IF (MATCHK (MXND, I1, I2, J1, J2, LXN)) THEN
         CONTINUE
C
C  TRY THE CURRENT I'S AND J'2 REVERSED
C
      ELSEIF (MATCHK (MXND, J1, J2, I1, I2, LXN)) THEN
         I1HOLD = I1
         I2HOLD = I2
         I1 = J1
         I2 = J2
         J1 = I1HOLD
         J2 = I2HOLD
         KOUNTL = NLOOP - KOUNTL - 2
      ELSE
C
C  TRY ONE STEP FORWARD AND BACKWARDS (NORMAL AND I'S AND J'S REVERSED)
C
         IFOR1 = I2
         IFOR2 = LNODES (3, I2)
         JFOR1 = J1
         JFOR2 = LNODES (2, J1)
C
         IBAC1 = LNODES (2, I1)
         IBAC2 = I1
         JBAC1 = LNODES (3, J2)
         JBAC2 = J2
C
         IF (MATCHK (MXND, IFOR1, IFOR2, JFOR2, JFOR1, LXN)) THEN
            I1 = IFOR1
            I2 = IFOR2
            J1 = JFOR2
            J2 = JFOR1
            KOUNTL = KOUNTL - 2
         ELSEIF (MATCHK (MXND, JFOR2, JFOR1, IFOR1, IFOR2, LXN)) THEN
            I1 = JFOR2
            I2 = JFOR1
            J1 = IFOR1
            J2 = IFOR2
            KOUNTL = NLOOP - KOUNTL
         ELSEIF (MATCHK (MXND, IBAC1, IBAC2, JBAC2, JBAC1, LXN)) THEN
            I1 = IBAC1
            I2 = IBAC2
            J1 = JBAC2
            J2 = JBAC1
            KOUNTL = KOUNTL + 2
         ELSEIF (MATCHK (MXND, JBAC2, JBAC1, IBAC1, IBAC2, LXN)) THEN
            I1 = JBAC2
            I2 = JBAC1
            J1 = IBAC1
            J2 = IBAC2
            KOUNTL = NLOOP - KOUNTL - 4
         ELSE
            ERR = .TRUE.
            GOTO 100
         ENDIF
C
      ENDIF
C
  100 CONTINUE
C
      RETURN
C
      END
