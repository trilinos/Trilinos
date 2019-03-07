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

C $Id: llist.f,v 1.1 1990/11/30 11:11:21 gdsjaar Exp $
C $Log: llist.f,v $
C Revision 1.1  1990/11/30 11:11:21  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]LLIST.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE LLIST (MS, ML, MAXNL, NS, NL, KNUM, LISTL, ILINE,
     &   ISIDE, NLPS, IFLINE, ILLIST, LCON, ISLIST, LINKS, LINKL, ERR)
C***********************************************************************
C
C  SUBROUTINE LLIST = PRODUCE LIST OF LINES FOR REGION
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     PERIM  = GENERATES PERIMETER OF THE REGION
C
C***********************************************************************
C
C     PRODUCE THE LIST OF  (PHYSICAL INDICES OF) LINES FOR  (PHYSICAL)
C     REGION KREG.  THIS LIST IS  (LISTL (I), I=1, NL).
C     *BACKWARDS* SIDES WILL BE REVERSED.
C     IT IS ASSUMED LINES ARE PROPERLY LISTED IN ORDER ON SIDE CARDS.
C     IF THEY ARE NOT,  PERIM WILL DIAGNOSE IT.
C     ERR = .TRUE. IF ERRORS WERE ENCOUNTERED.
C
C***********************************************************************
C
      DIMENSION ILINE (ML), LCON (3, ML)
      DIMENSION ISIDE (MS), NLPS (MS), IFLINE (MS), ILLIST (MS * 3)
      DIMENSION ISLIST (NS), LISTL (MAXNL)
      DIMENSION LINKL (2, ML), LINKS (2, MS)
C
      LOGICAL ERR, ADDLNK
C
      ERR = .TRUE.
      ADDLNK = .FALSE.
      IS = ISLIST (1)
C
C  FIRST SIDE
C
      IF (IS .EQ. 0) THEN
         RETURN
      ELSEIF (IS .LT. 0) THEN
         IFL1 = IABS (IS)
         ILL1 = IABS (IS)
         NEW = 1
         LISTL (NEW) = IFL1
      ELSE
         CALL LTSORT (MS, LINKS, IS, IPNTR, ADDLNK)
         I1 = IFLINE (IPNTR)
         I2 = I1 + NLPS (IPNTR) - 1
         IFL1 = ILLIST (I1)
         ILL1 = ILLIST (I2)
         NEW = 0
         DO 100 I = I1, I2
            NEW  =  NEW + 1
            LISTL (NEW) = ILLIST (I)
  100    CONTINUE
      ENDIF
      IF (NS .LE. 1) THEN
         NL = NEW
         ERR = .FALSE.
         RETURN
      ENDIF
C
C     SECOND SIDE
C
      IS2 = ISLIST (2)
      IF (IS2 .EQ. 0) THEN
         RETURN
      ELSEIF (IS2 .LT. 0) THEN
         IFL2 = IABS (IS2)
         ILL2 = IABS (IS2)
      ELSE
         CALL LTSORT (MS, LINKS, IS2, IPNTR, ADDLNK)
         I1 = IFLINE (IPNTR)
         I2 = I1 + NLPS (IPNTR) - 1
         IFL2 = ILLIST (I1)
         ILL2 = ILLIST (I2)
      ENDIF
C
C  DECIDE WHICH END OF SIDE ONE IS THE STARTING POINT
C
      CALL LTSORT (ML, LINKL, IFL2, IPNTR, ADDLNK)
      K1 = LCON (1, IPNTR)
      K2 = LCON (2, IPNTR)
      CALL LTSORT (ML, LINKL, ILL2, IPNTR, ADDLNK)
      K3 = LCON (1, IPNTR)
      K4 = LCON (2, IPNTR)
      CALL LTSORT (ML, LINKL, IFL1, IPNTR, ADDLNK)
      J1 = LCON (1, IPNTR)
      J2 = LCON (2, IPNTR)
      CALL LTSORT (ML, LINKL, ILL1, IPNTR, ADDLNK)
      J3 = LCON (1, IPNTR)
      J4 = LCON (2, IPNTR)
C
C  FIRST SIDE IN PROPER ORDER
C
      IF ( (J3 .EQ. K1) .OR. (J3 .EQ. K2) .OR. (J3 .EQ. K3)
     &   .OR. (J3 .EQ. K4) .OR. (J4 .EQ. K1) .OR. (J4 .EQ. K2)
     &   .OR. (J4 .EQ. K3) .OR. (J4 .EQ. K4)) THEN
         CONTINUE
C
C  FIRST SIDE NEEDS REVERSED
C
      ELSEIF ( (J1 .EQ. K1) .OR. (J1 .EQ. K2) .OR. (J1 .EQ. K3)
     &   .OR. (J1 .EQ. K4) .OR. (J2 .EQ. K1) .OR. (J2 .EQ. K2)
     &   .OR. (J2 .EQ. K3) .OR. (J2 .EQ. K4)) THEN
         CALL IREVER (LISTL, NEW)
C
C  CONNECTIVITY DOES NOT EXIST
C
      ELSE
         IF (IS2 .GT. 0) THEN
            CALL LTSORT (MS, LINKS, IS2, IPNTR, ADDLNK)
            WRITE ( * , 10000)KNUM, ISIDE (IPNTR)
         ELSE
            CALL LTSORT (ML, LINKL, IABS (IS2), IPNTR, ADDLNK)
            WRITE ( * , 10010)KNUM, ILINE (IPNTR)
         ENDIF
         RETURN
      ENDIF
C
      NL = NEW
      DO 120 KS = 2, NS
         I = LISTL (NL)
         CALL LTSORT (ML, LINKL, I, IPNTR, ADDLNK)
         J1 = LCON (1, IPNTR)
         J2 = LCON (2, IPNTR)
         IS = ISLIST (KS)
C
C     ADD NEW LINES TO LIST
C
         IF (IS .EQ. 0) THEN
            RETURN
         ELSEIF (IS .LT. 0) THEN
            IFL = IABS (IS)
            ILL = IABS (IS)
            NEW = NL + 1
            LISTL (NEW) = IABS (IS)
         ELSE
            CALL LTSORT (MS, LINKS, IS, IPNTR, ADDLNK)
            I1 = IFLINE (IPNTR)
            I2 = I1 + NLPS (IPNTR) - 1
            IFL = ILLIST (I1)
            ILL = ILLIST (I2)
            NEW = NL
            DO 110 I = I1, I2
               NEW = NEW + 1
               LISTL (NEW) = ILLIST (I)
  110       CONTINUE
         ENDIF
C
C     DETERMINE WHETHER THIS SIDE IS BACKWARDS
C
         CALL LTSORT (ML, LINKL, IFL, IPNTR, ADDLNK)
         K1 = LCON (1, IPNTR)
         K2 = LCON (2, IPNTR)
         CALL LTSORT (ML, LINKL, ILL, IPNTR, ADDLNK)
         K3 = LCON (1, IPNTR)
         K4 = LCON (2, IPNTR)
C
C  THIS SIDE IS IN PROPER ORDER
C
         IF ( (J1 .EQ. K1) .OR. (J1 .EQ. K2) .OR. (J2 .EQ. K1)
     &      .OR. (J2 .EQ. K2)) THEN
            CONTINUE
C
C  THIS SIDE NEEDS REVERSING
C
         ELSEIF ( (J1 .EQ. K3) .OR. (J1 .EQ. K4) .OR. (J2 .EQ. K3)
     &      .OR. (J2 .EQ. K4)) THEN
            CALL LTSORT (MS, LINKS, IS, IPNTR, ADDLNK)
            CALL IREVER (LISTL (NL + 1), NLPS (IPNTR))
C
C  CONNECTIVITY DOES NOT EXIST
C
         ELSE
            IF (IS .GT. 0) THEN
               CALL LTSORT (MS, LINKS, IS, IPNTR, ADDLNK)
               WRITE ( * , 10000)KNUM, ISIDE (IPNTR)
            ELSE
               CALL LTSORT (ML, LINKL, IABS (IS), IPNTR, ADDLNK)
               WRITE ( * , 10010)KNUM, ILINE (IPNTR)
            ENDIF
            RETURN
         ENDIF
C
         NL = NEW
C
  120 CONTINUE
C
C  SUCCESSFULL LINE LIST GENERATION
C
      ERR = .FALSE.
C
      RETURN
C
10000 FORMAT  (' IN REGION', I5, ', SIDE', I5, ' DOES NOT CONNECT TO',
     &   /, ' THE PREVIOUS SECTION OF THE PERIMETER')
10010 FORMAT  (' IN REGION', I5, ', LINE', I5, ' DOES NOT CONNECT TO',
     &   /, ' THE PREVIOUS SECTION OF THE PERIMETER')
C
      END
