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

C $Id: flmnmx.f,v 1.1 1990/11/30 11:07:38 gdsjaar Exp $
C $Log: flmnmx.f,v $
C Revision 1.1  1990/11/30 11:07:38  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]FLMNMX.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE FLMNMX (MXND, MLN, MAXPRM, LINKPR, KPERIM, LNODES,
     &   XN, YN, NLOOP, NODE, XMIN, XMAX, YMIN, YMAX, ERR)
C***********************************************************************
C
C  SUBROUTINE FLMNMX = SET MIN AND MAX FOR CURRENT FILL BOUNDARY
C
C***********************************************************************
C
      DIMENSION LNODES (MLN, MXND), XN (MXND), YN (MXND)
      DIMENSION LINKPR (3, MAXPRM)
C
      LOGICAL ERR
C
      KOUNT = 0
      INOW = NODE
      XMIN = XN (NODE)
      XMAX = XN (NODE)
      YMIN = YN (NODE)
      YMAX = YN (NODE)
C
  100 CONTINUE
C
      INOW = LNODES (3, INOW)
      IF (INOW .NE. NODE) THEN
C
         XMIN = MIN (XMIN, XN (INOW))
         YMIN = MIN (YMIN, YN (INOW))
         XMAX = MAX (XMAX, XN (INOW))
         YMAX = MAX (YMAX, YN (INOW))
C
         KOUNT = KOUNT + 1
C
         IF (KOUNT .GT. NLOOP) THEN
            CALL MESAGE('PROBLEMS IN FLMNMX WITH LOOP NOT CLOSING')
            ERR = .TRUE.
            GOTO 130
         ENDIF
         GOTO 100
      ENDIF
C
C  LOOP THROUGH ALL THE REMAINING PERIMETERS CHECKING FOR CROSSINGS
C
      IPERIM = KPERIM
  110 CONTINUE
      IPERIM = LINKPR (2, IPERIM)
      IF ((IPERIM .EQ. 0) .OR. (IPERIM .EQ. KPERIM)) GOTO 130
C
      KMAX = LINKPR (3, IPERIM)
      INOW = LINKPR (1, IPERIM)
      KOUNT = 0
C
  120 CONTINUE
      XMIN = MIN (XMIN, XN (INOW))
      YMIN = MIN (YMIN, YN (INOW))
      XMAX = MAX (XMAX, XN (INOW))
      YMAX = MAX (YMAX, YN (INOW))
      KOUNT = KOUNT + 1
      INOW = LNODES (3, INOW)
      IF (INOW .EQ. LINKPR (1, IPERIM)) GOTO 110
C
      IF (KOUNT. GT. KMAX + 1) THEN
         CALL MESAGE('PROBLEMS IN FLMNMX WITH LOOP NOT CLOSING')
         ERR = .TRUE.
         GOTO 130
      ENDIF
      GOTO 120
C
  130 CONTINUE
      RETURN
C
      END
