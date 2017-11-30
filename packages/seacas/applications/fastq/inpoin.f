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

C $Id: inpoin.f,v 1.1 1990/11/30 11:09:57 gdsjaar Exp $
C $Log: inpoin.f,v $
C Revision 1.1  1990/11/30 11:09:57  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]INPOIN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INPOIN (MP, N1, N18, JJ, X, Y, NHOLDP, IHOLDP, IPOINT,
     &   COOR, IPBOUN, LINKP, MERGE, NOROOM)
C***********************************************************************
C
C  SUBROUTINE INPOIN = ENTERS A POINT INTO THE DATABASE
C
C***********************************************************************
C
      DIMENSION IPOINT (MP), COOR (2, MP), IPBOUN (MP), LINKP (2, MP)
      DIMENSION IHOLDP (2, MP)
C
      LOGICAL NOROOM, MERGE, ADDLNK
C
      NOROOM = .TRUE.
      JHOLD = JJ
C
C  ZERO OUT THE LINK ARRAY IF NEEDED
C
      IF (JJ .GT. N18) THEN
         N18 = JJ
C
C  GET THE CORRECT NODE NUMBER IF MERGING
C
      ELSEIF (MERGE) THEN
         ADDLNK = .FALSE.
         CALL LTSORT (MP, LINKP, JJ, IPNTR, ADDLNK)
         IF (IPNTR .GT. 0) THEN
            IF (JHOLD .GT. NHOLDP)NHOLDP = JHOLD
            CALL LTSORT (MP, IHOLDP, JHOLD, IPNTR, ADDLNK)
            IF (IPNTR .GT. 0) THEN
               JJ = IPNTR
            ELSE
               JJ = N18 + 1
               N18 = N18 + 1
               WRITE ( * , 10000)JHOLD, JJ
               ADDLNK = .TRUE.
               CALL LTSORT (MP, IHOLDP, JHOLD, JJ, ADDLNK)
            ENDIF
         ENDIF
      ENDIF
C
C  INPUT THE POINT DATA
C
      N1 = N1 + 1
      J = N1
      IF (J .GT. MP)RETURN
      ADDLNK = .TRUE.
      CALL LTSORT (MP, LINKP, JJ, J, ADDLNK)
      IPOINT (J) = JJ
      COOR (1, J) = X
      COOR (2, J) = Y
      IPBOUN (J) = 0
      NOROOM = .FALSE.
      RETURN
C
10000 FORMAT ('  OLD POINT NO:', I5, '  TO NEW POINT NO:', I5)
      END
