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

C $Id: cclock.f,v 1.2 1991/03/21 15:44:23 gdsjaar Exp $
C $Log: cclock.f,v $
C Revision 1.2  1991/03/21 15:44:23  gdsjaar
C Changed all 3.14159... to atan2(0.0, -1.0)
C
c Revision 1.1.1.1  1990/11/30  11:04:21  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:04:19  gdsjaar
c Initial revision
c
C
CC* FILE: [.QMESH]CCLOCK.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CCLOCK (X, Y, N, CCW, ERR, INDETR)
C***********************************************************************
C
C  SUBROUTINE CCLOCK = DETERMINES IF THE PERIMETER OF A REGION IS STATED
C                      IN COUNTER-CLOCKWISE FASHION
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     PERIM = GENERATES THE PERIMETER OF A REGION
C
C***********************************************************************
C
C  VARIABLES USED:
C     CCW    = .TRUE. IF THE PERIMETER IS IN COUNTER-CLOCKWISE ORDER
C     ERR    = .TRUE. IF THE ORDER COULD NOT BE DETERMINED, OR IF AN
C              ERROR OCCURS CHECKING THE ORDER
C     N      = THE NUMBER OF NODES IN THE PERIMETER (MUST BE AT LEAST 3)
C
C***********************************************************************
C
      DIMENSION X (N), Y (N)
C
      LOGICAL CCW, ERR, INDETR
C
      ERR = .TRUE.
      INDETR = .FALSE.
      PI = ATAN2(0.0, -1.0)
      TWOPI = PI+PI
C
      IF (N .LT. 3) THEN
         CALL MESAGE ('PERIMETER MUST CONTAIN MORE THAN 3 NODES')
         GOTO 110
      ENDIF
C
      SPIRO = 0.0
      IF ( (X (1) .EQ. X (N)) .AND. (Y (1) .EQ. Y (N)) ) THEN
         CALL MESAGE ('PERIMETER CONTAINS DUPLICATE NODE LOCATIONS')
         GOTO 110
      ENDIF
      AGOLD = ATAN2 (Y (1) - Y (N), X (1) - X (N))
      DO 100 I = 1, N
         NEXT = I + 1
         IF (NEXT .GT. N)NEXT = 1
         IF ( (X (NEXT) .EQ. X (I)) .AND. (Y (NEXT) .EQ. Y (I)) ) THEN
            CALL MESAGE ('PERIMETER CONTAINS DUPLICATE NODE LOCATIONS')
            GOTO 110
         ENDIF
         AGNEW = ATAN2 (Y (NEXT) - Y (I), X (NEXT) - X (I))
         DIFF = AGNEW - AGOLD
         IF (DIFF .GT. PI) DIFF = DIFF-TWOPI
         IF (DIFF .LT. -PI) DIFF = DIFF+TWOPI
         IF (ABS (ABS (DIFF) - PI) .LT. 1.0E-3) THEN
            CALL MESAGE ('PERIMETER CONTAINS SWITCHBACKS')
            GOTO 110
         ENDIF
         SPIRO = SPIRO + DIFF
         AGOLD = AGNEW
  100 CONTINUE
      CCW = .TRUE.
      IF (SPIRO .LT. 0.0) CCW = .FALSE.
      IF ( (ABS (SPIRO) .LT. PI) .OR. (ABS (SPIRO) .GT. (3.*PI))) THEN
         INDETR = .TRUE.
         RETURN
      ENDIF
      ERR = .FALSE.
      RETURN
C
C  ERROR IN THE ROUTINE
C
  110 CONTINUE
      RETURN
C
      END
