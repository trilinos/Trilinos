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

C $Id: chkwt.f,v 1.1 1990/11/30 11:04:40 gdsjaar Exp $
C $Log: chkwt.f,v $
C Revision 1.1  1990/11/30 11:04:40  gdsjaar
C Initial revision
C
C
CC* FILE: [.RENUM]CHKWT.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CHKWT (MP, ML, MS, MLIST, MFLAG, NNUID, MXLPS, LINKP,
     &   LINKL, LINKS, NUID, IFLG, ILEN, IPTR, NLIST, IFLAG, LTYPE,
     &   LCON, NLPS, IFLINE, ILLIST, COOR, JPOINT, JSIDE, ILOC, JLOC,
     &   NIX, ILIST, XLIST, ADDLNK, ISPNT, ERR)
C***********************************************************************
C
C  SUBROUTINE CHKWT = CHECKS THE WEIGHTING COMBINATION FOR FLAGS TO
C                     MAKE SURE THE DEFINITION IS VALID.  IT RETURNS
C                     THE LOCATION OF THE BEGINNING POINT IN THE
C                     FLAG NODE LIST AS JLOC.  THE ARRAY ILIST AND
C                     XLIST ARE ALSO SET UP AS X VALUES OF THE SIDE
C                     AS IT MOVES FORWARD.
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     ADDWT = ADDS THE WEIGHTING FACTORS TO ANY NODES WITH
C             FLAGS CONTAINING WEIGHTS
C
C***********************************************************************
C
C
      DIMENSION IFLG (MFLAG), ILEN (MFLAG), IPTR (MFLAG), NLIST (MLIST)
      DIMENSION NUID (NNUID), LINKP (2, MP)
      DIMENSION LINKS (2, MS), LINKL (2, ML)
      DIMENSION COOR (2, MP), LTYPE (ML), LCON (3, ML)
      DIMENSION NLPS (MS), IFLINE (MS), ILLIST (MS*3)
      DIMENSION ILIST (MXLPS), XLIST (MXLPS)
C
      LOGICAL ERR, ADDLNK, ISPNT
C
      ERR = .TRUE.
C
C  CHECK TO MAKE SURE THE FLAG HAS NODES GENERATED ALONG IT
C
      DO 100 I = 1, MFLAG
         IF (IFLAG .EQ. IFLG (I)) THEN
            ILOC = I
            GOTO 110
         ENDIF
  100 CONTINUE
      WRITE (*, 10000)IFLAG
      RETURN
  110 CONTINUE
C
C  CHECK TO MAKE SURE THE BEGINNING POINT IS DEFINED
C
      CALL LTSORT (MP, LINKP, JPOINT, JPPNTR, ADDLNK)
      IF (JPPNTR .LE. 0) THEN
         WRITE (*, 10010)JPOINT, IFLAG
         RETURN
      ENDIF
C
C  CHECK TO MAKE SURE THE NODE AT THE BEGINNING POINT IS IN THE LIST
C  OF NODES FOR THE FLAG  -  RETURN THIS LIST LOCATION IF IT IS THERE
C
      JLOC = 0
      J1 = IPTR (ILOC)
      J2 = IPTR (ILOC) + ILEN (ILOC) - 1
      DO 120 I = J1, J2
         IF (NUID (NLIST (I)) .EQ. JPOINT) THEN
            JLOC = I
         ENDIF
  120 CONTINUE
      IF (JLOC .EQ. 0) THEN
         WRITE (*, 10020)JPOINT
         RETURN
      ENDIF
C
C  CHECK TO MAKE SURE THE SIDE OR LINE OR POINT IS DEFINED
C
      IF (ISPNT) THEN
         CALL LTSORT (MP, LINKP, JSIDE, IPNTR, ADDLNK)
         IF (IPNTR .LE. 0)WRITE (*, 10030)JSIDE, IFLAG
         RETURN
      ELSEIF (JSIDE .LT. 0) THEN
         CALL LTSORT (ML, LINKL, IABS (JSIDE), JLPNTR, ADDLNK)
         IF (JLPNTR .LE. 0) THEN
            WRITE (*, 10040)IABS (JSIDE), IFLAG
            RETURN
         ENDIF
      ELSE
         CALL LTSORT (MS, LINKS, JSIDE, JSPNTR, ADDLNK)
         IF (JSPNTR .LE. 0) THEN
            WRITE (*, 10050)JSIDE, IFLAG
            RETURN
         ENDIF
      ENDIF
C
C  CHECK TO MAKE SURE SIDE'S LINE DEFINITIONS ARE ALL THERE
C  AND ARE IN MONOTONICALLY INCREASING X ORDER
C
      IF (JSIDE .GT. 0) THEN
         J1 = IFLINE (JSPNTR)
         J2 = IFLINE (JSPNTR) + NLPS (JSPNTR) - 1
      ELSE
         J1 = 1
         J2 = 1
      ENDIF
      NIX = J2 - J1 + 2
      DO 130 J = J1, J2
         IF (JSIDE .GT. 0) THEN
            KK = ILLIST (J)
            CALL LTSORT (ML, LINKL, KK, IPNTR, ADDLNK)
            IF ((KK .LE. 0) .OR. (IPNTR .LE. 0)) THEN
               WRITE (*, 10060)JSIDE, KK
               RETURN
            ENDIF
         ELSE
            KK = IABS (JSIDE)
            IPNTR = JLPNTR
         ENDIF
         ILIST (J - J1 + 1) = IPNTR
C
C  CHECK TO MAKE SURE LINE'S POINT DEFINITIONS ARE ALL THERE
C
         IP1 = LCON (1, IPNTR)
         IP2 = LCON (2, IPNTR)
         IP3 = IABS (LCON (3, IPNTR))
         CALL LTSORT (MP, LINKP, IP1, IPNTR1, ADDLNK)
         CALL LTSORT (MP, LINKP, IP2, IPNTR2, ADDLNK)
         IF (IP3 .NE. 0) THEN
            CALL LTSORT (MP, LINKP, IABS (IP3), IPNTR3, ADDLNK)
         ELSE
            IP3 = 0
         ENDIF
C
         IF ((IP1 .LE. 0) .OR. (IPNTR1 .LE. 0)) THEN
            IF (JSIDE .GT. 0) THEN
               WRITE (*, 10070)KK, JSIDE, IP1
            ELSE
               WRITE (*, 10080)KK, IP1
            ENDIF
            RETURN
         ELSEIF ((IP2 .LE. 0) .OR. (IPNTR2 .LE. 0)) THEN
            IF (JSIDE .GT. 0) THEN
               WRITE (*, 10070)KK, JSIDE, IP2
            ELSE
               WRITE (*, 10080)KK, IP2
            ENDIF
            RETURN
         ELSEIF ((LTYPE (IPNTR) .NE. 1) .AND. ((IP3 .EQ. 0) .OR.
     &      (IPNTR3 .LE. 0))) THEN
            IF (JSIDE .GT. 0) THEN
               WRITE (*, 10070)KK, JSIDE, IP3
            ELSE
               WRITE (*, 10080)KK, IP3
            ENDIF
            RETURN
         ENDIF
  130 CONTINUE
C
C  GET THE XLIST VALUES SET UP FOR THIS SIDE
C
      IF (NIX .EQ. 2) THEN
         CALL LTSORT (MP, LINKP, LCON (1, ILIST (1)), IP1, ADDLNK)
         CALL LTSORT (MP, LINKP, LCON (2, ILIST (1)), IP2, ADDLNK)
         IF (COOR (1, IP1) .LT. COOR (1, IP2)) THEN
            XLIST (1) = COOR (1, IP1)
            XLIST (2) = COOR (1, IP2)
         ELSEIF (COOR (1, IP1) .GT. COOR (1, IP2)) THEN
            XLIST (1) = COOR (1, IP2)
            XLIST (2) = COOR (1, IP1)
         ELSE
            WRITE (*, 10090)JSIDE
            RETURN
         ENDIF
      ELSE
C
C  DEFINE VALUE OF KP,  THE BEGINNING CONNECTIVITY POINT
C
         K1 = LCON (1, ILIST (1))
         K2 = LCON (2, ILIST (1))
         L1 = LCON (1, ILIST (2))
         L2 = LCON (2, ILIST (2))
         KP = 0
         IF ((K1 .EQ. L1) .OR. (K1 .EQ. L2)) THEN
            KP = K1
            CALL LTSORT (MP, LINKP, K2, K2P, ADDLNK)
            XLIST (1) = COOR (1, K2P)
         ELSEIF ((K2 .EQ. L1) .OR. (K2 .EQ. L2)) THEN
            KP = K2
            CALL LTSORT (MP, LINKP, K1, K1P, ADDLNK)
            XLIST (1) = COOR (1, K1P)
         ELSE
            WRITE (*, 10100)ILLIST (J1), JSIDE
            RETURN
         ENDIF
C
C  NOWLOOP THROUGH THE REST OF THE LINES TO GET THE XLIST
C
         DO 140 M = 2, NIX - 1
            L1 = LCON (1, ILIST (M))
            L2 = LCON (2, ILIST (M))
            IF (KP .EQ. L1) THEN
               KP = L2
               CALL LTSORT (MP, LINKP, L1, L1P, ADDLNK)
               XLIST (M) = COOR (1, L1P)
            ELSEIF (KP .EQ. L2) THEN
               KP = L1
               CALL LTSORT (MP, LINKP, L2, L2P, ADDLNK)
               XLIST (M) = COOR (1, L2P)
            ELSE
               WRITE (*, 10100)ILLIST (J1 + M - 2), JSIDE
               RETURN
            ENDIF
  140    CONTINUE
         CALL LTSORT (MP, LINKP, KP, KPP, ADDLNK)
         XLIST (NIX) = COOR (1, KPP)
C
C  NOWCHECK TO MAKE SURE THE LINES ARE INCREASING MONOTONICALLY IN X
C
         DO 150 M = 2, NIX
            IF (XLIST (M) .LT. XLIST (M - 1)) THEN
               WRITE (*, 10090)JSIDE
               RETURN
            ENDIF
  150    CONTINUE
      ENDIF
      ERR = .FALSE.
      RETURN
C
10000 FORMAT (' FLAG NUMBER:', I5, ' HAS NOT BEEN USED IN GENERATING',
     &   /, ' THIS MESH - THUS NO WEIGHTING IS POSSIBLE')
10010 FORMAT (' THE BEGINNING POINT:', I5, ' FOR FLAG:', I5,
     &   ' IS NOT DEFINED ', /,
     &   ' NO BOUNDARY FLAG WEIGHTING WITH THIS POINT POSSIBLE')
10020 FORMAT (' THE BEGINNING POINT:', I5, ' IS NOT IN THE FLAG', /,
     &   ' NO BOUNDARY FLAG WEIGHTING POSSIBLE WITH THIS FLAG/POINT '//
     &   'SET.')
10030 FORMAT (' THE WEIGHTING POINT:', I5, ' FOR FLAG:', I5,
     &   ' IS NOT DEFINED', /,
     &   ' NO BOUNDARY FLAG WEIGHTING WITH THIS POINT POSSIBLE')
10040 FORMAT (' THE WEIGHTING LINE:', I5, ' FOR FLAG:', I5,
     &   ' IS NOT DEFINED', /,
     &   ' NO BOUNDARY FLAG WEIGHTING WITH THIS LINE POSSIBLE')
10050 FORMAT (' THE WEIGHTING SIDE:', I5, ' FOR FLAG:', I5,
     &   ' IS NOT DEFINED', /,
     &   ' NO BOUNDARY FLAG WEIGHTING WITH THIS SIDE POSSIBLE')
10060 FORMAT (' FOR SIDE:', I5, ' LINE:', I5, ' DOES NOT EXIST', /,
     &   ' NO BOUNDARY FLAG WEIGHTING WITH THIS SIDE POSSIBLE')
10070 FORMAT (' ON LINE:', I5, ' IN SIDE:', I5, ' POINT:', I5,
     &   ' DOES NOT EXIST', /,
     &   ' NO BOUNDARY FLAG WEIGHTING WITH THIS SIDE POSSIBLE')
10080 FORMAT (' ON LINE:', I5, ' POINT:', I5, ' DOES NOT EXIST', /,
     &   ' NO BOUNDARY FLAG WEIGHTING WITH THIS LINE POSSIBLE')
10090 FORMAT (' SIDE:', I5, ' IS NOT MONOTONICALLY INCREASING IN X', /,
     &   ' NO BOUNDARY FLAG WEIGHTING WITH THIS SIDE POSSIBLE')
10100 FORMAT (' LINE:', I5, ' DOES NOT CONNECT PROPERLY IN SIDE:', I5,
     &   /, ' NO BOUNDARY FLAG WEIGHTING WITH THIS SIDE POSSIBLE')
C
      END
