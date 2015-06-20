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

C $Id: resta.f,v 1.1 1990/11/30 11:14:47 gdsjaar Exp $
C $Log: resta.f,v $
C Revision 1.1  1990/11/30 11:14:47  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]RESTA.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE RESTA (MXND, XN, YN, NUID, LXK, KXL, NXL, LXN, KKK,
     &   KKKOLD, NAVAIL, IAVAIL, NNN, LIMIT, IREST, TILT, ERR, NOROOM)
C************************************************************************
C
C  SUBROUTINE RESTA  =  RESTRUCTURES THE MESH TO ELIMINATE WORST ELELMENTS
C
C***********************************************************************
C
C  NOTE:
C     A RECORD IS KEPT OF UP TO 25 OF THE CURRENT WORST CONDITION NUMBERS
C     AND THE WORST ELEMENT POSSIBLE IS RESTRUCTURED
C     UNTIL NO FURTHER RESTRUCTURING CAN BE DONE.
C
C***********************************************************************
C
      DIMENSION KCND(26), CND(26)
      DIMENSION NXL(2, 3*MXND), XN(MXND), YN(MXND), NUID(MXND)
      DIMENSION LXK(4, MXND), KXL(2, 3*MXND), LXN(4, MXND)
      DIMENSION NODES(4), ANGLES(4), SIDES(4)
C
      LOGICAL ERR, NOROOM, LSIDE, CCW, CAREA, DONE
C
      ERR = .FALSE.
C
C  CHECK FOR IMPENDING OVERFLOW
C
      IF (NAVAIL .LE. 1) THEN
         NOROOM = .TRUE.
         CALL MESAGE ('INSUFFICIENT STORAGE AVAILABLE IN RESTA')
         RETURN
      ENDIF
C
C  INITIALIZE
C
      NTAB = 0
      MAXTAB = 25
      CNDTOL = 2.0
      ASUM = 0.
      NSUM = 0
      CCW = .TRUE.
      CAREA = .FALSE.
      IREST = 0
C
      DO 110 K  =  KKKOLD  +  1, KKK
         IF (LXK(1, K) .GT. 0) THEN
            LSIDE = .FALSE.
C
C  GET THE ELEMENTS COND VALUE (BASED ON ANGLE AND SIDE LENGTH)
C
            CALL GNXKA (MXND, XN, YN, K, NODES, AREA, LXK, NXL, CCW)
            CALL QAAVAL (MXND, NODES, ANGLES, QRAT, DUMMY, XN, YN,
     &         CAREA)
            CALL CONDNO (MXND, NODES, QRAT, SRAT, COND, SIDES, XN, YN,
     &         LSIDE)
C
C  ADD UP THE NUMBER OF ANGLES < PI/2
C
            DO 100 I = 1, 4
               IF (ANGLES(I) .LE. 1.58) THEN
                  ASUM = ASUM + ANGLES(I)
                  NSUM = NSUM + 1
               ENDIF
  100       CONTINUE
C
C  ADD BAD ELEMENTS TO THE LIST
C
            IF (COND .GE. CNDTOL) THEN
               CND(NTAB + 1) = COND
               KCND(NTAB + 1) = K
               CALL BUBBLE (CND, KCND, NTAB, NTAB + 1)
               NTAB = MIN0(NTAB + 1, MAXTAB)
            ENDIF
         ENDIF
  110 CONTINUE
C
C  TILT IS THE AVERAGE VALUE IN DEGREES OF ANGLES < PI/2
C
      IF (NSUM .GT. 0) THEN
         TILT = (ASUM/FLOAT(NSUM))*57.2957795
      ELSE
         TILT = 90.
      ENDIF
      IF ((LIMIT .LE. 0) .OR. (NTAB .LE. 0)) RETURN
      CNDTOL = CND(NTAB)
C
C  TRY TO RESTRUCTURE ON THE 10 WORST ELEMENTS ONLY
C
  120 CONTINUE
      NTRY = MIN0(NTAB, 10)
      DO 130 IK = 1, NTRY
         IK1 = IK
         CALL RESTRY (MXND, KCND(IK), K2, LXK, NXL, KXL, LXN, XN, YN,
     &      NUID, NAVAIL, IAVAIL, NNN, DONE, ERR, NOROOM)
         IF (ERR) RETURN
         IF (DONE) GO TO 140
  130 CONTINUE
      RETURN
  140 CONTINUE
      IREST = IREST + 1
      IF (IREST .GE. LIMIT) RETURN
C
C  UPDATE THE TABLE (AFTER 1 RESTRUCTURE)
C
      CALL GNXKA (MXND, XN, YN, KCND(IK1), NODES, AREA, LXK, NXL, CCW)
      CALL QAAVAL (MXND, NODES, ANGLES, QRAT, DUMMY, XN, YN, CAREA)
      CALL CONDNO (MXND, NODES, QRAT, SRAT, COND1, SIDES, XN, YN, LSIDE)
      CND(IK1) = COND1
      DO 150 IK = 1, NTAB
         IK2 = IK
         IF (KCND(IK) .EQ. K2) GO TO 160
  150 CONTINUE
      IK2 = NTAB + 1
      NTAB = NTAB + 1
  160 CONTINUE
      CALL GNXKA (MXND, XN, YN, K2, NODES, AREA, LXK, NXL, CCW)
      CALL QAAVAL (MXND, NODES, ANGLES, QRAT, DUMMY, XN, YN, CAREA)
      CALL CONDNO (MXND, NODES, QRAT, SRAT, COND2, SIDES, XN, YN, LSIDE)
      CND(IK2) = COND2
      KCND(IK2) = K2
C
C  RE-SORT AND PRUNE
C
      CALL BUBBLE (CND, KCND, 1, NTAB)
      DO 170 I = 1, 2
         IF (CND(NTAB) .LT. CNDTOL)NTAB = NTAB - 1
  170 CONTINUE
      NTAB = MIN0(NTAB, MAXTAB)
      IF (NTAB .LE. 0) RETURN
      CNDTOL = CND(NTAB)
      GO TO 120
C
      END
