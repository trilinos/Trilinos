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

C $Id: trifix.f,v 1.3 1998/07/14 18:20:11 gdsjaar Exp $
C $Log: trifix.f,v $
C Revision 1.3  1998/07/14 18:20:11  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.2  1991/03/21 15:45:23  gdsjaar
C Changed all 3.14159... to atan2(0.0, -1.0)
C
c Revision 1.1.1.1  1990/11/30  11:17:19  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:17:17  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]TRIFIX.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE TRIFIX (MXND, MLN, XN, YN, ZN, NUID, LXK, KXL, NXL,
     &   LXN, NNN, LLL, KKK, NAVAIL, IAVAIL, ANGLE, LNODES, BNSIZE,
     &   NLOOP, DEV1, KREG, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAXZ, GRAPH,
     &   VIDEO, NOROOM, ERR)
C***********************************************************************
C
C  SUBROUTINE TRIFIX = CHECKS ALL ELEMENTS FOR ANY TRIANGULAR SHAPED
C                      LONG ELEMENT AND DELETES THEM WHEN
C                      FOUND AND POSSIBLE
C
C***********************************************************************
C
      DIMENSION ANGLE (MXND), BNSIZE (2, MXND), LNODES (MLN, MXND)
      DIMENSION NODES(4)
      DIMENSION LXK(4, MXND), NXL(2, 3*MXND), KXL(2, 3*MXND)
      DIMENSION LXN(4, MXND), XN(MXND), YN(MXND), ZN(MXND), NUID(MXND)
C
      CHARACTER*3 DEV1
      LOGICAL ERR, DONE, GRAPH, REDO, CCW
      LOGICAL VIDEO, NOROOM
C
      PI = ATAN2(0.0, -1.0)
      TWOPI = 2.0 * PI
C
      ERR = .FALSE.
      DONE = .FALSE.
      CCW = .TRUE.
      KMAX = 30
      KOUNT = 0
C
C  TOLERANCE IS SET AT 150 DEGREES
C
      TOLER = 2.6179939
C
  100 CONTINUE
      KOUNT = KOUNT + 1
      IF (KOUNT .GT. KMAX) GOTO 140
      REDO = .FALSE.
C
      DO 130 KELEM = 1, KKK
         IF (LXK (1, KELEM) .GT. 0) THEN
            CALL GNXKA (MXND, XN, YN, KELEM, NODES, AREA, LXK, NXL, CCW)
            DO 110 I = 1, 4
               I1 = NODES (I)
               IF (I .EQ. 1) THEN
                  I0 = NODES (4)
                  I2 = NODES (2)
               ELSEIF (I .EQ. 4) THEN
                  I0 = NODES (3)
                  I2 = NODES (1)
               ELSE
                  I0 = NODES (I - 1)
                  I2 = NODES (I + 1)
               ENDIF
C
               ANG1 = ATAN2 (YN (I0) - YN (I1), XN (I0) - XN (I1))
               IF (ANG1 .LT. 0.) ANG1 = ANG1 + TWOPI
               ANG2 = ATAN2 (YN (I2) - YN (I1), XN (I2) - XN (I1))
               IF (ANG2 .LT. 0.) ANG2 = ANG2 + TWOPI
               ANG = ANG1 - ANG2
               IF (ANG .LT. 0.) ANG = ANG + TWOPI
C
               CALL LONGEL (MXND, MLN, LNODES, XN, YN, NUID, LXK, KXL,
     &            NXL, LXN, NNN, NAVAIL, IAVAIL, I1, KELEM, ANG, TOLER,
     &            I0, I2, KREG, XMIN, XMAX, YMIN, YMAX, KKK, LLL,
     &            DONE, GRAPH, VIDEO, NOROOM, ERR, KKKADD)
               IF ((NOROOM) .OR. (ERR)) GOTO 140
C
               IF (DONE) THEN
C
                  IF ((GRAPH) .AND. (.NOT. VIDEO)) THEN
                     DIST = MAX (ABS(XN (I0) - XN (I1)),
     &                  ABS(XN (I2) - XN (I1)), ABS(YN (I0) - YN (I1)),
     &                  ABS(YN (I2) - YN (I1))) * 3.
                     XMIN = XN (I1) - DIST
                     XMAX = XN (I1) + DIST
                     YMIN = YN (I1) - DIST
                     YMAX = YN (I1) + DIST
                  ENDIF
                  IF (VIDEO) THEN
                     CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX,
     &                  YMIN, YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
                     CALL SNAPIT (3)
                  ENDIF
C
                  CALL FILSMO (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL,
     &               LXN, LLL, NNN, NNN, LNODES, BNSIZE, NLOOP, XMIN,
     &               XMAX, YMIN, YMAX, ZMIN, ZMAX, DEV1, KREG)
                  IF ((GRAPH) .OR. (VIDEO)) THEN
                     CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX,
     &                  YMIN, YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
                     IF (VIDEO) CALL SNAPIT (3)
                  ENDIF
                  DONE = .FALSE.
                  REDO = .TRUE.
                  GOTO 120
               ENDIF
C
  110       CONTINUE
  120       CONTINUE
         ENDIF
  130 CONTINUE
C
      IF (REDO) GOTO 100
  140 CONTINUE
C
      RETURN
C
      END
