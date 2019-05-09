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

C $Id: gnxka.f,v 1.2 2000/11/13 15:41:35 gdsjaar Exp $
      SUBROUTINE GNXKA (MXND, XN, YN, K, NODES, AREA, LXK, NXL, CCW)
C***********************************************************************
C  SUBROUTINE GNXKA = GENERATES A LIST OF THE FOUR NODES ASSOCIATED WITH
C                     ELEMENT K
C
C***********************************************************************
C  VARIABLES USED:
C     CCW    = .TRUE. IF LIST IS TO BE IN CCW ORDER AND AREA DEFINED
C    (Changed to always put in order and calculate area)
C
C***********************************************************************
C
      REAL XN (MXND), YN (MXND)
      INTEGER NODES(4)
      INTEGER LXK(4, MXND), NXL(2, 3 * MXND)

      LOGICAL CCW

      AREA = 0.0
      DO 10 I = 1, 4
        NODES (I) = 0
 10   CONTINUE

C... Let line 1 be the base line
      L = LXK(1, K)

C... Null List
      IF (L .LE. 0) THEN
        RETURN
      ENDIF

      NODES(1) = NXL(1, L)
      NODES(2) = NXL(2, L)

C... Find other ends of the two sides

      DO 110 I = 2, 4
         L  = LXK(I, K)
         M1 = NXL(1, L)
         M2 = NXL(2, L)

         IF (M1 .EQ. NODES(1)) THEN
           NODES(4) = M2
         ELSE IF (M2 .EQ. NODES(1)) THEN
           NODES(4) = M1
         END IF

         IF (M1 .EQ. NODES(2)) THEN
           NODES(3) = M2
         ELSE IF (M2 .EQ. NODES(2)) THEN
           NODES(3) = M1
         END IF

  110 CONTINUE

C... Compute signed area
      AREA = 0.5 *
     *  ((XN(NODES(3)) - XN(NODES(1))) *
     *   (YN(NODES(4)) - YN(NODES(2))) -
     &   (YN(NODES(3)) - YN(NODES(1))) *
     *   (XN(NODES(4)) - XN(NODES(2))))

      IF (AREA .LT. 0.) THEN
C ... Clockwise case  -  reverse the order
        NTMP = NODES(2)
        NODES(2) = NODES(4)
        NODES(4) = NTMP
        AREA = -AREA
      ENDIF

      RETURN
      END
