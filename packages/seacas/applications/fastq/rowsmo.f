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

C $Id: rowsmo.f,v 1.5 2004/01/22 14:25:22 gdsjaar Exp $
C $Log: rowsmo.f,v $
C Revision 1.5  2004/01/22 14:25:22  gdsjaar
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
C Revision 1.4  1998/11/24 20:45:09  gdsjaar
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
C Revision 1.3  1998/07/14 18:19:59  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.2  1991/03/21 15:45:16  gdsjaar
C Changed all 3.14159... to atan2(0.0, -1.0)
C
c Revision 1.1.1.1  1990/11/30  11:15:12  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:15:10  gdsjaar
c Initial revision
c
C
CC* FILE: [.PAVING]ROWSMO.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE ROWSMO (MXND, MLN, XN, YN, ZN, LXK, KXL, NXL, LXN, NNN,
     &   WFAC, WFAC2, NIT, EPS, RO, NNN2, LNODES, BNSIZE, LLL, GRAPH,
     &   XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, DEV1, KREG)
C***********************************************************************
C
C  SUBROUTINE ROWSMO = SMOOTHS AN ADDED ROW DURING FILLING USING THE
C                      ISOPARAMETRIC SMOOTHER WITH SPECIAL CONSIDERATION
C                      GIVEN TO THE 2-LINE NODES (ROW CORNERS)
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
      DIMENSION AREA(20)
C
      DIMENSION KLIST(20), NODES(4)
      DIMENSION XN(MXND), YN(MXND), ZN(MXND)
      DIMENSION LXN(4, MXND), NXL(2, 3*MXND)
      DIMENSION LXK(4, MXND), KXL(2, 3*MXND)
      DIMENSION LINES(20), LNODES (MLN, MXND), BNSIZE (2, MXND)
C
      LOGICAL BIG, CCW, ERR, GRAPH, NEAR2L, TEST, AVER
C
      CHARACTER*3 DEV1
C
      PI = ATAN2(0.0, -1.0)

      nodes(1) = -1
      nodes(2) = -1
      nodes(3) = -1
      nodes(4) = -1

      IF (RO .LT. 0.01) RO = 1.
      DRO = 1.0
      VRO = 1.0
      EPS2 = EPS * RO
      TEST = .FALSE.
      AVER = .TRUE.
C
C  ITERATION LOOP
C
      DO 160 IT = 1, NIT
         IF (IT .EQ. NIT) THEN
            CALL MESAGE ('THE ROWSMO ROUTINE IS NOT CONVERGING')
         ENDIF
         BIG = .FALSE.
C
C  NODE LOOP
C
         NBEGIN = MAX0 (NNN2 - 1, 1)
         NEND = NNN + 1
         DO 150 J = NBEGIN, NEND
            IF (J .EQ. NEND) THEN
               NODE = LNODES (3, NNN)
            ELSEIF (J .EQ. NBEGIN) THEN
               NODE = LNODES (2, NBEGIN)
            ELSE
               NODE = J
            ENDIF
C
C  SKIP CONTINUATION LINES, EXTERIOR BOUNDARY LINES, AND NODES
C  THAT ARE ON THE INTERIOR
C
            IF (NODE .gt. 0 .AND.
     &         (LXN(1, NODE) .GT. 0) .AND. (LXN(2, NODE) .GT. 0) .AND.
     &         (LNODES (4, NODE) .EQ. - 1) ) THEN
C
C  FIND ELEMENTS AND LINES ATTACHED TO NODE
C
               CALL GKXN (MXND, KXL, LXN, NODE, KS, KLIST, ERR)
               CALL GETLXN (MXND, LXN, NODE, LINES, NL, ERR)
C
               SUMX = 0.0
               SUMY = 0.0
C
C  PERFORM AN AREA PULL AND LAPLACIAN
C  ON ANY NODE ATTACHED TO A 2-LINE NODE
C
C               TWOL = .FALSE.
C               NEAR2L = .FALSE.
C               IF (LXN (3, NODE) . EQ. 0) THEN
C                  TWOL = .TRUE.
C                  NFROM = 0
C               ELSE
C                  IF (NL .EQ. 3) THEN
C                     DO 100 IL = 1, NL
C                        ILL = LINES (IL)
C                        IF (NXL (1, ILL) .EQ. NODE) THEN
C                           NTEST = NXL (2, ILL)
C                        ELSEIF (NXL (2, ILL) .EQ. NODE) THEN
C                           NTEST = NXL (1, ILL)
C                        ELSE
C                           CALL MESAGE ('** PROBLEMS IN ROWSMO **')
C                           GOTO 110
C                        ENDIF
C                        NODES(IL) = NTEST
C                        IF (LXN (3, NTEST) .EQ. 0) THEN
C
C  MAKE SURE THAT THE OTHER END OF THE 2-LINE NODE HAS ONLY 3 LINES
C
C                           IF (LXN (1, NTEST) .EQ. ILL) THEN
C                              LTEST = LXN (2, NTEST)
C                           ELSEIF (LXN (2, NTEST) .EQ. ILL) THEN
C                              LTEST = LXN (1, NTEST)
C                           ELSE
C                              CALL MESAGE ('** PROBLEMS IN ROWSMO **')
C                              GOTO 110
C                           ENDIF
C                           IF (NXL (1, LTEST) .EQ. NTEST) THEN
C                              NTEST = NXL (2, LTEST)
C                           ELSEIF (NXL (2, LTEST) .EQ. NTEST) THEN
C                              NTEST = NXL (1, LTEST)
C                           ELSE
C                              CALL MESAGE ('** PROBLEMS IN ROWSMO **')
C                              GOTO 110
C                           ENDIF
C                           IF ((LXN (3, NTEST) .GT. 0) .AND.
C     &                        (LXN (4, NTEST) .EQ. 0)) THEN
C                              NEAR2L = .TRUE.
C                           ENDIF
C                        ENDIF
C  100                CONTINUE
C  110                CONTINUE
C                     NFROM = NODES(2)
C                  ELSE
C                     NFROM = 0
C                  ENDIF
C               ENDIF
C
               NEAR2L = .FALSE.
               IF ((NEAR2L) .OR. (TEST)) THEN
                  THETA1 = ATAN2 (YN (NODES (3)) - YN (NODES (1)),
     &               XN (NODES (3)) - XN (NODES (1)) ) + PI / 2.0
                  THETA2 = ATAN2 (YN (NODES (3)) - YN (NODES (2)),
     &               XN (NODES (3)) - XN (NODES (2)) ) + PI / 2.0
                  DET =  - COS (THETA1) * SIN (THETA2) + COS (THETA2) *
     &               SIN (THETA1)
                  X11 = 0.5 * (XN (NODES (1)) + XN (NODES (3)))
                  Y11 = 0.5 * (YN (NODES (1)) + YN (NODES (3)))
                  X21 = 0.5 * (XN (NODES (2)) + XN (NODES (3)))
                  Y21 = 0.5 * (YN (NODES (2)) + YN (NODES (3)))
                  R =  (- SIN (THETA2) * (X21 - X11) + COS (THETA2) *
     &               (Y21 - Y11)) / DET
                  XNEW = X11 + R * COS (THETA1)
                  YNEW = Y11 + R * SIN (THETA1)
                  XDEL = XNEW - XN (NODE)
                  YDEL = YNEW - YN (NODE)
C
C  PERFORM AN ISOPARAMETRIC SMOOTH ON OTHER NODES
C
               ELSE
                  DO 120 KL = 1, KS
                     CCW = .FALSE.
                     KK = KLIST(KL)
                     CALL GNXKA (MXND, XN, YN, KK, NODES, AREA(KL), LXK,
     &                  NXL, CCW)
C
                     DO 100 IN = 1, 4
                        IF (NODES(IN) .EQ. NODE) THEN
                           J1 = IN + 1
                           IF (J1 .GT. 4) J1 = 1
                           GO TO 110
                        END IF
  100                CONTINUE
  110                CONTINUE
                     J2 = J1 + 1
                     IF (J2 .GT. 4) J2 = 1
                     J3 = J2 + 1
                     IF (J3 .GT. 4) J3 = 1
C
                     SUMX = SUMX + XN(NODES(J1)) + XN(NODES(J3))
     &                  - WFAC * XN(NODES(J2))
                     SUMY = SUMY + YN(NODES(J1)) + YN(NODES(J3))
     &                  - WFAC * YN(NODES(J2))
  120             CONTINUE
                  SUMX = SUMX/(FLOAT(KS) * (2.0 - WFAC))
                  SUMY = SUMY/(FLOAT(KS) * (2.0 - WFAC))
                  XDEL = (RO * ( SUMX - XN (NODE) ))
                  YDEL = (RO * ( SUMY - YN (NODE) ))
C
                  CALL GETFRM (MXND, LINES, NL, NXL, NODE,
     &               LNODES (2, NODE), LNODES (3, NODE), NFROM)
                  IF (NFROM .GT. 0)
     &               THEN
C
C  FACTOR IN THE LENGTH CONSTANT (GENERATED LENGTH) OF THE NODE
C
                     DIST0 = BNSIZE (1,NODE) * BNSIZE (2,NODE)
                     XDIST = XDEL + XN (NODE) - XN (NFROM)
                     YDIST = YDEL + YN (NODE) - YN (NFROM)
                     DIST1 = SQRT (XDIST **2 + YDIST **2)
                     DFACT = (DIST0 / DIST1) * DRO
                     SUMX = XN (NFROM) + XDIST * DFACT
                     SUMY = YN (NFROM) + YDIST * DFACT
                     XDEL = SUMX - XN (NODE)
                     YDEL = SUMY - YN (NODE)
C
C  FACTOR IN THE EQUAL ANGLE VECTORS
C
                     IF (LNODES (2, NODE) .NE. LNODES (3, NODE)) THEN
                        CALL EQLANG (MXND, XN, YN, LXN, NODE,
     &                     LNODES (2, NODE), LNODES (3, NODE), NFROM,
     &                     DIST0, VRO, VX, VY)
                        IF (AVER) THEN
                           XDEL = (XDEL + VX) * .5
                           YDEL = (YDEL + VY) * .5
                        ELSE
                           XDEL = VX
                           YDEL = VY
                        ENDIF
                     ENDIF
C
                  ENDIF
C
               ENDIF
C
C  NOW CHECK THAT THE ROW IS NOT BENDING OVER ON ITSELF WITH THIS SMOOTH
C
               IF (LXN (4, NODE) .EQ. 0) CALL INVERT_FQ (MXND, MLN, XN,
     &            YN, ZN, LXK, KXL, NXL, LXN, LLL, LNODES, XMIN, XMAX,
     &            YMIN, YMAX, ZMIN, ZMAX, DEV1, KREG, NODE, XDEL, YDEL)
C
C  REDEFINE THIS NODE'S COORDINATES
C  AND PLOT THE NEW NODE AND LINES
C
               IF ((XDEL * XDEL + YDEL * YDEL) .GT. EPS2) BIG = .TRUE.
               IF (GRAPH) THEN
                  CALL LCOLOR ('BLACK')
                  DO 130 II = 1, NL
                     IDRAW = LINES(II)
                     NODE1 = NXL (1, IDRAW)
                     NODE2 = NXL (2, IDRAW)
                     CALL D2NODE (MXND, XN, YN, NODE1, NODE2)
  130             CONTINUE
                  CALL LCOLOR ('WHITE')
                  CALL SFLUSH
               ENDIF
C
               XN(NODE) = XN(NODE) + XDEL
               YN(NODE) = YN(NODE) + YDEL
C
               IF (GRAPH) THEN
                  DO 140 II = 1, NL
                     IDRAW = LINES(II)
                     NODE1 = NXL (1, IDRAW)
                     NODE2 = NXL (2, IDRAW)
                     CALL D2NODE (MXND, XN, YN, NODE1, NODE2)
  140             CONTINUE
                  CALL SFLUSH
               ENDIF
C
C  CHECK FOR CONVERGENCE
C
            ENDIF
  150    CONTINUE
C
C  IF NO SIGNIFICANT MOVEMENTS OCCURRED,  RETURN
C
         IF (.NOT.BIG) RETURN
  160 CONTINUE
C
C  NOW SMOOTH THE INTERIOR
C
      RETURN
      END
