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

C $Id: filsmo.f,v 1.2 1998/07/14 18:18:57 gdsjaar Exp $
C $Log: filsmo.f,v $
C Revision 1.2  1998/07/14 18:18:57  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:07:25  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:07:23  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]FILSMO.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE FILSMO (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN,
     &   LLL, NNN, NNN2, LNODES, BNSIZE, NLOOP, XMIN, XMAX, YMIN, YMAX,
     &   ZMIN, ZMAX, DEV1, KREG)
C***********************************************************************
C
C  SUBROUTINE FILSMO = MESH SMOOTHING DONE BY ISOPARAMETRIC/EQUAL
C                      ANGULAR SMOOTHING OF THE ADDED INTERIOR (FREE)
C                      BOUNDARY ROW AND THEN A LENGTH-WEIGHTED/EQUAL
C                      ANGULAR BOUNDARY LAPLACIAN OF THE INTERIOR NODES.
C                      THE FREE BOUNDARY IS FINALLY SMOOTHED AGAIN.
C
C***********************************************************************
C
C  VARIABLES USED:
C     WFAC = WEIGTH (0. = LAPLACIAN, 1. = ISOPARAMETRIC)
C     NIT  = THE MAX NUMBER OF ITERATIONS TO DO.
C     EPS  = MINIMUM DISTANCE NODES MUST MOVE TO CONTINUE ITERATIONS
C     RO   = AN UNDER- OR OVER-RELAXATION FACTOR (NORMALLY 1.0)
C
C***********************************************************************
C
      COMMON /TIMING/ TIMEA, TIMEP, TIMEC, TIMEPC, TIMEAJ, TIMES
C
      DIMENSION XN(MXND), YN(MXND), ZN(MXND)
      DIMENSION LXN(4, MXND), NXL(2, 3*MXND)
      DIMENSION LXK(4, MXND), KXL(2, 3*MXND)
      DIMENSION LINES(20), LNODES (MLN, MXND), BNSIZE (2, MXND)
C
      LOGICAL BIG, ERR, GRAPH, DONE
C
      CHARACTER*3 DEV1
C
      CALL GETIME (TIME1)
      GRAPH = .FALSE.
      DONE = .FALSE.
      WT = 10.
C
      NIT = MAX0 (5 * NLOOP, 40)
      TOL = .03
      VRO = 1.
      RO = 1.
      WFAC = 1.0
      WFAC2 = .5
      CALL MNORM  (MXND,  XN,  YN,  NXL,  LLL,  STDLEN)
      EPS  =  TOL * STDLEN
      IF (RO .LT. 0.01) RO = 1.
      EPS2 = (EPS * RO)**2
C
C  FIRST SMOOTH THE ADDED ROW
C
      IF (NLOOP .GT. 0) THEN
         CALL ROWSMO (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN, NNN,
     &      WFAC, WFAC2, NIT, EPS, RO, NNN2, LNODES, BNSIZE, LLL,
     &      GRAPH, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, DEV1, KREG)
      ENDIF
C
C  NOW SMOOTH THE INTERIOR NODES
C
C  ITERATION LOOP
C
      DO 140 IT = 1, NIT
         BIG = .FALSE.
C
C  NODE LOOP
C
         DO 130 NODE = 1, NNN
            IF ( (LXN (1, NODE) .GT. 0) .AND.
     &         (LXN (2, NODE) .GT. 0) .AND.
     &         (LNODES (4, NODE) .EQ. - 2) ) THEN
               DONE = .TRUE.
               FX = 0.
               FY = 0.
               SL = 0.
               VL = 0.
C
C  LOOP THRU ALL LINES CONNECTED TO NODE
C
               CALL GETLXN (MXND, LXN, NODE, LINES, KOUNT, ERR)
               IF (ERR) GOTO 150
               DO 100 IL = 1, KOUNT
                  L = LINES (IL)
                  NEND = NXL (1, L) + NXL (2, L) - NODE
                  DX = XN (NEND) - XN (NODE)
                  DY = YN (NEND) - YN (NODE)
                  AL = SQRT (DX * DX + DY * DY)
C
C  CHECK FOR A BOUNDARY NODE AT THE OTHER END
C  OF THE LINE - TRY TO AVERAGE ANGULAR ERRORS WITH THE BOUNDARY WHERE
C  POSSIBLE - THIS MEANS ADDING IN AN EXTRA VECTOR TO PULL THE NODE
C  BACK TO WHERE IT OUGHT TO BE TO BE AT EQUAL ANGLES
C
                  IF (LXN (2, NEND) .LT. 0) THEN
                     CALL SETN02 (MXND, NXL, LXK, KXL, L, NEND, NODE,
     &                  N0, N2)
                     CALL EQLANG (MXND, XN, YN, LXN, NODE, N0, N2,
     &                  NEND, AL, VRO, VXDEL, VYDEL)
                     VL = SQRT (VXDEL * VXDEL + VYDEL * VYDEL)
                     FX = FX + (VXDEL * WT * VL)
                     FY = FY + (VYDEL * WT * VL)
                     SL = SL + VL * WT
                  ENDIF
                  FX = FX + DX * AL
                  FY = FY + DY * AL
                  SL = SL + AL
  100          CONTINUE
C
C  MOVE THE NODE
C
               DELX = RO * FX/SL
               DELY = RO * FY/SL
C
C  ERASE THE NODE'S LINES IF GRAPH IS ON
C
               IF (GRAPH) THEN
                  CALL LCOLOR('BLACK')
                  DO 110 II = 1, KOUNT
                     IDRAW = LINES(II)
                     NODE1 = NXL (1, IDRAW)
                     NODE2 = NXL (2, IDRAW)
                     CALL D2NODE (MXND, XN, YN, NODE1, NODE2)
  110             CONTINUE
                  CALL LCOLOR ('WHITE')
               ENDIF
C
               XN (NODE) = XN (NODE)+DELX
               YN (NODE) = YN (NODE)+DELY
C
C  REPLOT THE NODE'S LINES IF GRAPH IS ON
C
               IF (GRAPH) THEN
                  DO 120 II = 1, KOUNT
                     IDRAW = LINES(II)
                     NODE1 = NXL (1, IDRAW)
                     NODE2 = NXL (2, IDRAW)
                     CALL D2NODE (MXND, XN, YN, NODE1, NODE2)
  120             CONTINUE
                  CALL SFLUSH
               ENDIF
               IF (DELX ** 2 + DELY ** 2 .GT. EPS2) BIG = .TRUE.
            ENDIF
  130    CONTINUE
         IF (.NOT.BIG) GOTO 150
  140 CONTINUE
  150 CONTINUE
C
C  NOW RESMOOTH THE ADDED ROW IF THE MESH HAS CHANGED INTERNALLY
C
      IF ((NLOOP .GT. 0) .AND. (DONE)) THEN
         CALL ROWSMO (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN, NNN,
     &      WFAC, WFAC2, NIT, EPS, RO, NNN2, LNODES, BNSIZE, LLL,
     &      GRAPH, XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, DEV1, KREG)
      ENDIF
C
C  NOW RESET ALL THE NODES AS BEING SMOOTHED
C
      DO 160 I = 1, NNN
         LNODES (4, I) = IABS (LNODES (4, I))
  160 CONTINUE
C
      CALL GETIME (TIME2)
      TIMES = TIMES + TIME2 - TIME1
      RETURN
C
      END
