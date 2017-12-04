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

C $Id: getang.f,v 1.5 2004/01/26 17:28:18 gdsjaar Exp $
C $Log: getang.f,v $
C Revision 1.5  2004/01/26 17:28:18  gdsjaar
C Removed several unused variables from getang subroutine.
C
C Initialized a variable
C
C Revision 1.4  2004/01/23 21:05:26  gdsjaar
C Removed integer*8 statement incorrectly checked in from AMD port
C
C Revision 1.3  2004/01/22 14:25:22  gdsjaar
C Attempt to fix strange problem on x86_64 AMD Opteron system using
C Portland Group 5.1-3 compilers. The getang function would work
C correctly if compiled with no optimization and in debug mode, but
C would crash if compiled optimized. The location of the crash was not
C in a place that made any sense that something was wrong.
C
C After much trial and error, it was found that adding a 'SAVE'
C statement at the beginning of the file fixed the problem.
C
C Also cleaned out some unused parameters being passed to the function.
C
C Revision 1.2  1991/03/21 15:44:47  gdsjaar
C Changed all 3.14159... to atan2(0.0, -1.0)
C
c Revision 1.1.1.1  1990/11/30  11:07:54  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:07:53  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]GETANG.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GETANG (MXND, MLN, XN, YN, LNODES, LXK, KXL, NXL,
     &  LXN, I, J, K, ANGLE, ERR)
C***********************************************************************
C
C SUBROUTINE GETANG = RETURNS THE CCW ANGLE FROM A VECTOR DRAWN
C                     FROM NODE J TO K TO A VECTOR DRAWN
C                      FROM NODE J TO I
C
C***********************************************************************
C
      DIMENSION LNODES (MLN, *)
      DIMENSION XN (*), YN (*)
      DIMENSION LXN(4, *), NXL(2, *)
      DIMENSION LXK(4, *), KXL(2, *)
C
      LOGICAL CORNP, SIDEP, DISCTP, ERR
C

C ... The save statement was added during debugging on the AMD Opteron
C     system using the pgf77 5.1-3 compiler.  Without the save, the
C     code coredumps at line 120 if optimized....  Could not track
C     down a problem, but adding the SAVE did work...
      
      SAVE

      TWOPI = 2.0 * ATAN2(0.0, -1.0)

      IOPP = 0
      KOPP = 0
      I1 = 0
      KK1 = 0
C
C  CHECK FOR NODES ON TOP OF EACH OTHER
C
      IF (((XN (J) .EQ. XN (K)) .AND. (YN (J) .EQ. YN (K)) ) .OR.
     &  ( (XN (I) .EQ. XN (J)) .AND. (YN (I) .EQ. YN (J)) ) .OR.
     &  ( (XN (I) .EQ. XN (K)) .AND. (YN (I) .EQ. YN (K)) ) ) THEN
        ANGLE = 0.
        GOTO 220
      ENDIF
C
      V1 = ATAN2 (YN (K)-YN (J), XN (K)-XN (J))
      IF (V1 .LT. 0.) V1 = V1 + TWOPI
      V2 = ATAN2 (YN (I)-YN (J), XN (I)-XN (J))
      IF (V2 .LT. 0.) V2 = V2 + TWOPI
      ANGLE = V2 - V1
      IF (ANGLE .LT. 0.) ANGLE = ANGLE + TWOPI
C
C  NOW CHECK TO MAKE SURE THAT THE ANGLE HAS NOT CROSSED THE PREVIOUS
C  ELEMENTS SIDES
C
      L1 = LNODES (5, I)
      L2 = LNODES (5, J)
      K1 = KXL (1, L1)
      K2 = KXL (1, L2)
      IF (K1 .EQ. K2) GOTO 210
C
C  SEE IF L2 CROSSES INTO K1 - FIRST GET THE NODE OPPOSITE I
C  AND THEN CHECK THE ANGLE FROM VECTOR J TO K AND VECTOR
C  J TO IOPP AGAINST THE INTERNAL ANGLE - SMALLER AND IT HAS
C CROSSED OVER.
C
      IF (K1 .NE. 0) THEN
        DO 100 II = 1, 4
          LTEST = LXK (II, K1)
          IF (LTEST .NE. L1) THEN
            IF (NXL (1, LTEST) .EQ. J) THEN
              IOPP = NXL (2, LTEST)
              L3 = LTEST
              GOTO 110
            ELSEIF (NXL (2, LTEST) .EQ. J) THEN
              IOPP = NXL (1, LTEST)
              L3 = LTEST
              GOTO 110
            ENDIF
          ENDIF
 100    CONTINUE
C
        CALL MESAGE ('** PROBLEMS IN GETANG GETTING IOPP **')
        ERR = .TRUE.
        GOTO 220
 110    CONTINUE

C
C  NOW TEST FOR CROSS-OVER
C
        V2OPP = ATAN2 (YN (IOPP) - YN (J), XN (IOPP) - XN (J))
        IF (V2OPP .LT. 0.) V2OPP = V2OPP + TWOPI
        ANGLE1 = V2OPP - V2
        IF (ANGLE1 .LT. 0.) ANGLE1 = ANGLE1 + TWOPI
        ANGLE2 = V2OPP - V1
        IF (ANGLE2 .LT. 0.) ANGLE2 = ANGLE2 + TWOPI
        IF (ANGLE2 .LE. ANGLE1) THEN
          ANGLE = ANGLE - TWOPI
          GOTO 210
        ENDIF
      END IF
C
C  SEE IF L2 CROSSES INTO K1 - FIRST GET THE NODE OPPOSITE K
C  AND THEN CHECK THE ANGLE FROM VECTOR J TO K AND VECTOR
C  J TO IOPP AGAINST THE INTERNAL ANGLE - SMALLER AND IT HAS
C  CROSSED OVER.
C
 120  CONTINUE
      IF (K2 .EQ. 0) GOTO 210
      DO 130 II = 1, 4
        LTEST = LXK (II, K2)
        IF (LTEST .NE. L2) THEN
          IF (NXL (1, LTEST) .EQ. J) THEN
            KOPP = NXL (2, LTEST)
            GOTO 140
          ELSEIF (NXL (2, LTEST) .EQ. J) THEN
            KOPP = NXL (1, LTEST)
            GOTO 140
          ENDIF
        ENDIF
 130  CONTINUE
C
      CALL MESAGE ('** PROBLEMS IN GETANG GETTING KOPP **')
      ERR = .TRUE.
      GOTO 220
 140  CONTINUE
C
C  NOW TEST FOR CROSS-OVER
C
      V1OPP = ATAN2 (YN (KOPP) - YN (J), XN (KOPP) - XN (J))
      IF (V1OPP .LT. 0.) V1OPP = V1OPP + TWOPI
      ANGLE1 = V1 - V1OPP
      IF (ANGLE1 .LT. 0.) ANGLE1 = ANGLE1 + TWOPI
      ANGLE2 = V2 - V1OPP
      IF (ANGLE2 .LT. 0.) ANGLE2 = ANGLE2 + TWOPI
      IF (ANGLE2 .LE. ANGLE1) THEN
        ANGLE = ANGLE - TWOPI
        GOTO 210
      ENDIF
C
C  NOW CHECK TO MAKE SURE THAT THE ANGLE HAS NOT CROSSED OVER TWO
C  ELEMENT SIDES IF THE NODE IS ATTACHED TO 5 OR MORE LINES
C
      IF (LXN (4, J) .LT. 0) THEN
        K3 = KXL (1, L3) + KXL (2, L3) - K1
        IF (K3 .EQ. K2) GOTO 210
C
C  SEE IF L2 CROSSES INTO K3 - FIRST GET THE NODE OPPOSITE J
C  AND THEN CHECK THE ANGLE FROM VECTOR J TO K AND VECTOR
C  J TO IOPP AGAINST THE INTERNAL ANGLE - SMALLER AND IT HAS
C CROSSED OVER.
C
        IF (K3 .EQ. 0) GOTO 120
        DO 150 II = 1, 4
          LTEST = LXK (II, K3)
          IF (LTEST .NE. L3) THEN
            IF (NXL (1, LTEST) .EQ. J) THEN
              IOPP3 = NXL (2, LTEST)
              GOTO 160
            ELSEIF (NXL (2, LTEST) .EQ. J) THEN
              IOPP3 = NXL (1, LTEST)
              GOTO 160
            ENDIF
          ENDIF
 150    CONTINUE
C
        CALL MESAGE ('** PROBLEMS IN GETANG GETTING IOPP3 **')
        ERR = .TRUE.
        GOTO 220
 160    CONTINUE

C
C  NOW TEST FOR CROSS-OVER
C
        V3OPP = ATAN2 (YN (IOPP3) - YN (J), XN (IOPP3) - XN (J))
        IF (V3OPP .LT. 0.) V3OPP = V3OPP + TWOPI
        ANGLE1 = V3OPP - V2
        IF (ANGLE1 .LT. 0.) ANGLE1 = ANGLE1 + TWOPI
        ANGLE2 = V3OPP - V1
        IF (ANGLE2 .LT. 0.) ANGLE2 = ANGLE2 + TWOPI
        IF (ANGLE2 .LE. ANGLE1) THEN
          ANGLE = ANGLE - TWOPI
          GOTO 210
        ENDIF
      ENDIF
C
C  NOW CHECK FOR AN INVERTED THREE NODE ANGLE - VERY SPECIAL
C  CASE THAT FALLS THROUGH THE PREVIOUS CHECK
C
      IF (KOPP .EQ. IOPP) THEN
        DO 170 II = 1, 4
          LTEST = LXK (II, K1)
          IF ( (NXL (1, LTEST) .EQ. IOPP) .AND.
     &      (NXL (2, LTEST) .NE. J) ) THEN
            I1 = NXL (2, LTEST)
            GOTO 180
          ELSEIF ( (NXL (2, LTEST) .EQ. IOPP) .AND.
     &        (NXL (1, LTEST) .NE. J) ) THEN
            I1 = NXL (1, LTEST)
            GOTO 180
          ENDIF
 170    CONTINUE
C
        CALL MESAGE ('** PROBLEMS IN GETANG GETTING I1 **')
        GOTO 220
 180    CONTINUE
        DO 190 II = 1, 4
          LTEST = LXK (II, K2)
          IF ( (NXL (1, LTEST) .EQ. KOPP) .AND.
     &      (NXL (2, LTEST) .NE. J) ) THEN
            KK1 = NXL (2, LTEST)
            GOTO 200
          ELSEIF ( (NXL (2, LTEST) .EQ. KOPP) .AND.
     &        (NXL (1, LTEST) .NE. J) ) THEN
            KK1 = NXL (1, LTEST)
            GOTO 200
          ENDIF
 190    CONTINUE
C
        CALL MESAGE ('** PROBLEMS IN GETANG KK1 **')
        GOTO 220
 200    CONTINUE

C
C  NOW TEST FOR INVERSION
C
        VVJ = ATAN2 (YN (J) - YN (KOPP), XN (J) - XN (KOPP))
        IF (VVJ .LT. 0.) VVJ = VVJ + TWOPI
        VVI1 = ATAN2 (YN (I1) - YN (KOPP), XN (I1) - XN (KOPP))
        IF (VVI1 .LT. 0.) VVI1 = VVI1 + TWOPI
        VVK1 = ATAN2 (YN (KK1) - YN (KOPP), XN (KK1) - XN (KOPP))
        IF (VVK1 .LT. 0.) VVK1 = VVK1 + TWOPI
        ANGLE1 = VVI1 - VVK1
        IF (ANGLE1 .LT. 0.) ANGLE1 = ANGLE1 + TWOPI
        ANGLE2 = VVJ - VVK1
        IF (ANGLE2 .LT. 0.) ANGLE2 = ANGLE2 + TWOPI
        IF (ANGLE2 .GT. ANGLE1) THEN
          ANGLE = ANGLE - TWOPI
        ENDIF
      ENDIF
C
C  GET THE RIGHT CLASSIFICATION
C
 210  CONTINUE
      IF (CORNP (ANGLE)) THEN
        IF (SIDEP (ANGLE)) THEN
          LNODES (6, J) = 2
        ELSE
          LNODES (6, J) = 1
        ENDIF
      ELSEIF (SIDEP (ANGLE)) THEN
        IF (DISCTP (ANGLE)) THEN
          LNODES (6, J) = 4
        ELSE
          LNODES (6, J) = 3
        ENDIF
      ELSE
        LNODES (6, J) = 5
      ENDIF
C
 220  CONTINUE
      RETURN
C
      END
