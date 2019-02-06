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

C $Id: nps.f,v 1.1 1990/11/30 11:12:52 gdsjaar Exp $
C $Log: nps.f,v $
C Revision 1.1  1990/11/30 11:12:52  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]NPS.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE NPS (ML, MS, MNNPS, NS, ISLIST, NINT, IFLINE, NLPS,
     &   ILLIST, LINKL, LINKS, NNPS, ERR)
C***********************************************************************
C
C  SUBROUTINE NPS = GIVES A LIST OF THE NUMBER OF PERIMETER NODES
C                   ON EACH OF A REGION'S SIDES
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     QMESH = GENERATES QUAD ELEMENTS
C
C***********************************************************************
C
C  VARIABLES USED:
C     NNPS  = ARRAY OF NUMBER OF NODES PER SIDE
C     CCW   = .TRUE. IF THE SIDE IS ORIENTED CCW
C     KS    = COUNTER OF THE NUMBER OF SIDES
C
C***********************************************************************
C
      DIMENSION NNPS (MNNPS), ISLIST (NS), LINKL (2, ML), LINKS (2, MS)
      DIMENSION NLPS (MS), IFLINE (MS), ILLIST (MS * 3), NINT (ML)
C
      LOGICAL ERR, ADDLNK
C
      ERR = .TRUE.
      ADDLNK = .FALSE.
C
      KS = 0
      DO 110 I = 1, NS
         IF (ISLIST (I) .LT. 0) THEN
            KS = KS + 1
            IL = IABS (ISLIST (I))
            CALL LTSORT (ML, LINKL, IL, IPNTR, ADDLNK)
            NNPS (KS) = IABS (NINT (IPNTR)) + 1
         ELSEIF (ISLIST (I) .GT. 0) THEN
            CALL LTSORT (MS, LINKS, ISLIST (I), IPNTR, ADDLNK)
            J1 = IFLINE (IPNTR)
            J2 = J1 + NLPS (IPNTR) - 1
            KS = KS + 1
            NNPS (KS) = 0
            DO 100 J = J1, J2
               IL = ILLIST (J)
               CALL LTSORT (ML, LINKL, IL, IPNTR, ADDLNK)
               NNPS (KS) = NNPS (KS) + IABS (NINT (IPNTR))
  100       CONTINUE
            NNPS (KS) = NNPS (KS) + 1
         ELSE
            RETURN
         ENDIF
  110 CONTINUE
      ERR = .FALSE.
C
      RETURN
C
      END
