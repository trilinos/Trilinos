C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Government retains certain rights in this software.
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

C $Id: inside.f,v 1.1 1990/11/30 11:10:14 gdsjaar Exp $
C $Log: inside.f,v $
C Revision 1.1  1990/11/30 11:10:14  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]INSIDE.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INSIDE (MS, N3, N4, N20, JJ, IIN, IFOUND, ISIDE, NLPS,
     &   IFLINE, ILLIST, LINKS, NHOLDS, IHOLDS, MERGE, NOROOM)
C***********************************************************************
C
C  SUBROUTINE INSIDE = INPUTS A SIDE INTO THE DATABASE
C
C***********************************************************************
C
      DIMENSION ISIDE (MS), NLPS (MS), IFLINE (MS), ILLIST (3 * MS)
      DIMENSION LINKS (2, MS)
      DIMENSION IIN (IFOUND), IHOLDS (2, MS)
C
      LOGICAL MERGE, NOROOM, ADDLNK
C
      IZ = 0
      NOROOM = .TRUE.
C
C  ZERO OUT THE LINK ARRAY IF NEEDED
C
      IF (JJ .GT. N20) THEN
         N20 = JJ
C
C  FIND THE CORRECT LINE NUMBER IF MERGING
C
      ELSEIF (MERGE) THEN
         ADDLNK = .FALSE.
         CALL LTSORT (MS, LINKS, JJ, IPNTR, ADDLNK)
         JHOLD = JJ
         IF (IPNTR .GT. 0) THEN
            IF (JHOLD .GT. NHOLDS)NHOLDS = JHOLD
            CALL LTSORT (MS, IHOLDS, JHOLD, IPNTR, ADDLNK)
            IF (IPNTR .GT. 0) THEN
               JJ = IPNTR
            ELSE
               JJ = N20 + 1
               N20 = N20 + 1
               WRITE ( * , 10000)JHOLD, JJ
               ADDLNK = .TRUE.
               CALL LTSORT (MS, IHOLDS, JHOLD, JJ, ADDLNK)
            ENDIF
         ENDIF
      ENDIF
C
C  INPUT THE SIDE DATA INTO THE DATABASE
C
      N3 = N3 + 1
      J = N3
      IF (J .GT. MS)RETURN
      ADDLNK = .TRUE.
      CALL LTSORT (MS, LINKS, JJ, J, ADDLNK)
      ISIDE (J) = JJ
      IFLINE (J) = N4 + 1
      DO 100 I = 1, IFOUND
         JJ = IIN (I)
         IF (JJ .EQ. 0)GOTO 110
         N4 = N4 + 1
         IF (N4 .GT. MS * 3)RETURN
         ILLIST (N4) = JJ
  100 CONTINUE
  110 CONTINUE
      NLPS (J) = N4 - IFLINE (J) + 1
      IF (NLPS (J) .LT. 1) THEN
         WRITE ( * , 10010)J
         CALL LTSORT (MS, LINKS, ISIDE (J), IZ, ADDLNK)
      ENDIF
C
      NOROOM = .FALSE.
      RETURN
C
10000 FORMAT ('   OLD SIDE NO:', I5, '   TO NEW SIDE NO:', I5)
10010 FORMAT (' SIDE:', I5, ' HAS LESS THAN ONE LINE', / ,
     &   ' THIS SIDE WILL NOT BE INPUT INTO DATABASE')
      END
