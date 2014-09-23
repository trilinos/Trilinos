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

C $Id: pcross.f,v 1.4 2000/11/13 15:39:05 gdsjaar Exp $
C $Log: pcross.f,v $
C Revision 1.4  2000/11/13 15:39:05  gdsjaar
C Cleaned up unused variables and labels.
C
C Removed some real to int conversion warnings.
C
C Revision 1.3  1999/06/17 19:02:22  gdsjaar
C Fixed several problems related to holes.  In several places, a
C nonpositive integer was being used to index into an array.  This seems
C to fix all of those cases.  I'm not sure if I fixed the true cause of
C these errors or just the symptom though...
C
C Revision 1.2  1998/07/14 18:19:30  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:13:11  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:13:10  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]PCROSS.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE PCROSS (MXND, MXCORN, MLN, MXLOOP, MAXPRM, NUID,
     &   XN, YN, ZN, LXK, KXL, NXL, LXN, ANGLE, LNODES, BNSIZE, LINKPR,
     &   KPERIM, NODE, NODE1, NODE2, KKKOLD, LLLOLD, NNNOLD, IAVAIL,
     &   NAVAIL, DONE, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, DEV1, LLL,
     &   KKK, NNN, LCORN, NCORN, NLOOP, NEXTN1, KLOOP, GRAPH, VIDEO,
     &   KREG, NOROOM, ERR)
C***********************************************************************
C
C  SUBROUTINE PCROSS = CHECKS TO SEE IF ANY PERIMETERS CROSS AND HOOKS
C                      THEM TOGETHER IF THEY DO
C
C***********************************************************************
C
      COMMON /TIMING/ TIMEA, TIMEP, TIMEC, TIMEPC, TIMEAJ, TIMES
C
      DIMENSION XN (MXND), YN (MXND), ZN (MXND), NUID (MXND)
      DIMENSION LXK (4, MXND), KXL (2, 3*MXND)
      DIMENSION NXL (2, 3*MXND), LXN (4, MXND)
      DIMENSION ANGLE (MXND), LNODES (MLN, MXND), BNSIZE (2, MXND)
      DIMENSION LCORN (MXCORN)
      DIMENSION NLOOP (MXLOOP), NEXTN1 (MXLOOP), LINKPR (3, MAXPRM)
C
      CHARACTER*3 DEV1
C
      LOGICAL DONE,  ERR, NOROOM, DONE1
      LOGICAL GRAPH, BOK, LCROSS, LMATCH
      LOGICAL VIDEO, PMATCH
C
C  FIND THE FIRST OVERLAPPING LINE STARTING AT THE CURRENT NODE
C
      CALL GETIME (TIME1)
      ERR = .FALSE.
      PMATCH = .TRUE.
C
  100 CONTINUE
      if (node1 .eq. 0) return
      N1 = NODE1
      KOUNT = 0
C
  110 CONTINUE
      N0 = LNODES (2, N1)
      N2 = LNODES (3, N1)
      N3 = LNODES (3, N2)
      KOUNT = KOUNT + 1
C
C  CHECK FOR COMPLETION
C
      IF ((N1 .EQ. NODE2) .AND. (KOUNT .GT. 1)) THEN
         GOTO 140
      ELSEIF (KOUNT .GT. NLOOP (1) + 1) THEN
         CALL MESAGE ('** PROLEMS WITH LOOP CLOSING IN PCROSS **')
         ERR = .TRUE.
         GOTO 140
      ENDIF
C
C  LOOP THROUGH ALL THE REMAINING PERIMETERS CHECKING FOR CROSSINGS
C
      IPERIM = LINKPR (2, KPERIM)
  120 CONTINUE
C
      IF (IPERIM .EQ. KPERIM) THEN
         N1 = N2
         GOTO 110
      ENDIF
C
      KOUNT2 = 0
      N1TEST = LINKPR (1, IPERIM)
C
  130 CONTINUE
      N0TEST = LNODES (2, N1TEST)
      N2TEST = LNODES (3, N1TEST)
      N3TEST = LNODES (3, N2TEST)
C
      CALL INTSCT (XN(N1), YN(N1), XN(N2), YN(N2), XN(N1TEST),
     &   YN(N1TEST), XN(N2TEST), YN(N2TEST), U, W, LCROSS)
      IF (.NOT. LCROSS) THEN
         N1TEST = N2TEST
         KOUNT2 = KOUNT2 + 1
         IF (N1TEST .EQ. LINKPR (1, IPERIM)) THEN
            IPERIM = LINKPR (2, IPERIM)
            GOTO 120
         ELSEIF (KOUNT2 .GT. LINKPR (3, IPERIM)) THEN
            CALL MESAGE ('** PROBLEMS IN PCROSS WITH UNCLOSED '//
     &         'PERIMETER **')
            ERR = .TRUE.
            GOTO 140
         ENDIF
         GOTO 130
      ENDIF
C
C  AN INTERSECTION HAS OCCURRED.
C  GET THE BEST SEAM FROM THIS INTERSECTION
C
      IF ((GRAPH) .OR. (VIDEO)) THEN
         CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX,
     &      YMIN, YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
         IF (VIDEO) CALL SNAPIT (2)
         IF (GRAPH) THEN
            CALL LCOLOR ('YELOW')
            CALL D2NODE (MXND, XN, YN, N1, N2)
            CALL D2NODE (MXND, XN, YN, N1TEST, N2TEST)
            CALL LCOLOR ('WHITE')
            CALL SFLUSH
         ENDIF
      ENDIF
      CALL MATCH2 (MXND, MLN, XN, YN, NXL, LXN, LNODES, ANGLE, N0, N1,
     &   N2, N3, N0TEST, N1TEST, N2TEST, N3TEST, I1, I2, J1, J2,
     &   KOUNTL, LMATCH, LINKPR (3, IPERIM), NODE, U, W, NLOOP (1),
     &   PMATCH, ERR)
      IF (ERR) GOTO 140
      IF (GRAPH) THEN
         CALL LCOLOR ('PINK ')
         CALL D2NODE (MXND, XN, YN, I1, I2)
         CALL D2NODE (MXND, XN, YN, J1, J2)
         CALL LCOLOR ('WHITE')
         CALL SFLUSH
      ENDIF
      IF (.NOT. LMATCH) THEN
         IF (N1TEST .EQ. LINKPR (1, IPERIM)) THEN
            IPERIM = LINKPR (2, IPERIM)
            GOTO 120
         ELSEIF (KOUNT2 .GT. LINKPR (3, IPERIM)) THEN
            CALL MESAGE ('** PROBLEMS IN PCROSS WITH UNCLOSED '//
     &         'PERIMETER **')
            ERR = .TRUE.
            GOTO 140
         ENDIF
         GOTO 130
      ENDIF
C
C  NOW CHECK TO SEE IF THE ATTACHMENT WOULD CAUSE
C  LINES ON THE BOUNDARY TO CROSS
C
      CALL BCROSS (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN, LNODES,
     &   I1, I2, J1, J2, NLOOP(1), BOK, LLL, XMIN, XMAX, YMIN, YMAX,
     &   ZMIN, ZMAX, DEV1, KREG, ERR)
      IF (ERR) GOTO 140
      IF (.NOT. BOK) THEN
         N1TEST = N2TEST
         KOUNT2 = KOUNT2 + 1
         IF (N1TEST .EQ. LINKPR (1, IPERIM)) THEN
            IPERIM = LINKPR (2, IPERIM)
            GOTO 120
         ELSEIF (KOUNT2 .GT. LINKPR (3, IPERIM)) THEN
            CALL MESAGE ('** PROBLEMS IN PCROSS WITH UNCLOSED '//
     &         'PERIMETER **')
            ERR = .TRUE.
            GOTO 140
         ENDIF
         GOTO 130
      ENDIF
C
C  NOW THAT THE APPROPRIATE COLLAPSE HAS BEEN FOUND, THE TWO LINES
C  MUST BE JOINED AND THE PERIMETER LINKS RESTABLISHED
C
      CALL SEW2 (MXND, MLN, NUID, LXK, KXL, NXL, LXN, LNODES,
     &   IAVAIL, NAVAIL, LLL, KKK, NNN, I1, I2, J1, J2, NOROOM, ERR)
      IF ((NOROOM) .OR. (ERR)) GOTO 140
C
C  UPDATE THE CURRENT NODE
C
      IF (J1 .EQ. NODE) THEN
         NDUM = NODE
         NODE = I2
         IF (NODE1 .EQ. NDUM) NODE1 = I2
         IF (NODE2 .EQ. NDUM) NODE2 = I2
      ELSEIF (J2 .EQ. NODE) THEN
         NDUM = NODE
         NODE = I1
         IF (NODE1 .EQ. NDUM) NODE1 = I1
         IF (NODE2 .EQ. NDUM) NODE2 = I1
      ENDIF
C
      NLOOP (1) = NLOOP (1) + LINKPR (3, IPERIM) - 2
      LINKPR (3, KPERIM) = NLOOP (1)
      JPERIM = LINKPR (2, IPERIM)
      IF (JPERIM .EQ. KPERIM) KPERIM = IPERIM
      LINKPR (1, IPERIM) = LINKPR (1, JPERIM)
      LINKPR (2, IPERIM) = LINKPR (2, JPERIM)
      LINKPR (3, IPERIM) = LINKPR (3, JPERIM)
      IF (LINKPR (2, KPERIM) .EQ. KPERIM) LINKPR (2, KPERIM) = 0
C
C  NOW SMOOTH AND PLOT THE CURRENT MESH
C
      NNN2 = 1
      CALL GETIME (TIME2)
      TIMEPC = TIMEPC + TIME2 - TIME1
      CALL FILSMO  (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN, LLL, NNN,
     &   NNN2, LNODES, BNSIZE, NLOOP (1), XMIN, XMAX, YMIN, YMAX, ZMIN,
     &   ZMAX, DEV1, KREG)
      CALL GETIME (TIME1)
      IF ((GRAPH) .OR. (VIDEO)) THEN
         CALL RPLOTL (MXND, XN, YN, ZN, NXL, XMIN, XMAX,
     &      YMIN, YMAX, ZMIN, ZMAX, LLL, DEV1, KREG)
         IF (VIDEO) CALL SNAPIT (2)
         IF (GRAPH) THEN
            CALL LCOLOR ('YELOW')
            CALL D2NODE (MXND, XN, YN, I1, I2)
            CALL LCOLOR ('WHITE')
            CALL SFLUSH
         ENDIF
      ENDIF
C
C  NOW TRY TO PINCH THE CONNECTION
C
      CALL LUPANG (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN, NLOOP (1),
     &   ANGLE, LNODES, I2, LLL, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX,
     &   DEV1, KREG, ERR)
      IF (ERR) GOTO 140

      CALL GETIME (TIME2)
      TIMEPC = TIMEPC + TIME2 - TIME1
      CALL PINCH (MXND, MXCORN, MLN, NUID, XN, YN, ZN, LXK, KXL, NXL,
     &   LXN, ANGLE, LNODES, BNSIZE, NODE, NLOOP (1), KKKOLD, LLLOLD,
     &   NNNOLD, IAVAIL, NAVAIL, DONE1, XMIN, XMAX, YMIN, YMAX, ZMIN,
     &   ZMAX, DEV1, LLL, KKK, NNN, LCORN, NCORN, NODE1, NODE2, GRAPH,
     &   VIDEO, KREG, NOROOM, ERR)
      IF ((NOROOM) .OR. (ERR)) GOTO 140
      IF (LINKPR(2, KPERIM) .NE. 0) GO TO 100
      CALL GETIME (TIME1)
C
  140 CONTINUE
      CALL GETIME (TIME2)
      TIMEPC = TIMEPC + TIME2 - TIME1
      RETURN
C
      END
