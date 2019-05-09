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

C $Id: qaaval.f,v 1.2 1991/03/21 15:45:11 gdsjaar Exp $
C $Log: qaaval.f,v $
C Revision 1.2  1991/03/21 15:45:11  gdsjaar
C Changed all 3.14159... to atan2(0.0, -1.0)
C
c Revision 1.1.1.1  1990/11/30  11:14:03  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:14:02  gdsjaar
c Initial revision
c
C
CC* FILE: [.QMESH]QAAVAL.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE QAAVAL (MXND, NODES, ANGLES, QRAT, AREA, XN, YN, CAREA)
C***********************************************************************
C
C  SUBROUTINE QAAVAL = CALCULATES THE INTERIOR ANGLES OF A QUAD AND
C                      THE RATIO OF LARGEST TO SMALLEST ANGLE
C
C***********************************************************************
C
C  VARIABLES USED:
C     NODES  = THE FOUR NODES OF THE ELEMENT IN CCW ORDER
C     ANGLES = THE FOUR INTERIOR ANGLES IN THE ORDER 4-1, 1-2, 2-3, 3-4
C     QRAT   = THE RATIO OF THE LARGEST TO THE SMALLEST ANGLE
C     CAREA  = .TRUE. IF THE AREA IS TO BE CALCULATED
C     AREA   = ELEMENT'S AREA
C
C***********************************************************************
C
      DIMENSION NODES (4), ANGLES (4), AG (4), XN (MXND), YN (MXND)
C
      LOGICAL CAREA
C
      PI = ATAN2(0.0, -1.0)
      TWOPI = PI + PI
C
      IF (CAREA) THEN
         N1 = NODES (1)
         N2 = NODES (2)
         N3 = NODES (3)
         N4 = NODES (4)
         AREA = 0.5 *  ( (XN (N3) - XN (N1)) *  (YN (N4) - YN (N2))
     &      -  (YN (N3) - YN (N1)) *  (XN (N4) - XN (N2)))
      ENDIF
C
      DO 100 I = 1, 4
         IF (I .EQ. 4) THEN
            J = 1
         ELSE
            J = I + 1
         ENDIF
         N1 = NODES (I)
         N2 = NODES (J)
         IF ( (XN (N2) .EQ. XN (N1)) .AND. (YN (N2) .EQ. YN (N1))) THEN
            QRAT = 1.E + 10
            RETURN
         ELSE
            AG (J) = ATAN2 (YN (N2) - YN (N1), XN (N2) - XN (N1))
         ENDIF
  100 CONTINUE
C
      DO 110 J = 1, 4
         IF (J .EQ. 1) THEN
            I = 4
         ELSE
            I = J - 1
         ENDIF
         DIFF = AG (J) - AG (I)
         IF (DIFF .GT. PI) THEN
            DIFF = DIFF - TWOPI
         ELSEIF (DIFF .LT.  - PI) THEN
            DIFF = DIFF + TWOPI
         ENDIF
         ANGLES (J) = PI - DIFF
  110 CONTINUE
C
      QMIN = ANGLES (1)
      QMAX = ANGLES (1)
      DO 120 I = 2, 4
         QMIN = AMIN1 (QMIN, ANGLES (I))
         QMAX = AMAX1 (QMAX, ANGLES (I))
  120 CONTINUE
      IF (QMIN .GT. 0.) THEN
         QRAT = QMAX / QMIN
      ELSE
         QRAT = 1.0E10
      ENDIF
C
      RETURN
C
      END
