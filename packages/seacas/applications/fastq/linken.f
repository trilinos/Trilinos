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

C $Id: linken.f,v 1.2 2000/11/13 15:39:04 gdsjaar Exp $
C $Log: linken.f,v $
C Revision 1.2  2000/11/13 15:39:04  gdsjaar
C Cleaned up unused variables and labels.
C
C Removed some real to int conversion warnings.
C
C Revision 1.1.1.1  1990/11/30 11:11:10  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:11:09  gdsjaar
c Initial revision
c
C
CC* FILE: [.MAIN]LINKEN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE LINKEN (MDIM, JJ, IFLAG1, IFLAG, IFLIST, NEPS, LIST,
     &   LINKF, LINKE, IBOUN, ADDLNK)
C***********************************************************************
C
C  SUBROUTINE LINKEN = LINKS ENTITIES IN BOUNDARY FLAG LISTS
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     LINKBC = LINKS UP ALL BOUNDARY FLAG LISTS
C
C***********************************************************************
C
C  VARIABLES USED:
C     IFLAG  = THE ARRAY OF FLAGS
C     IFLIST = THE FIRST ENTITY IN LIST TO BE ASSOCIATED WITH A FLAG
C     NEPS   = THE NUMBER OF ENTITIES IN LIST THAT GO WITH A FLAG
C     LIST   = THE LIST OF ENTITIES
C     LINK   = THE LINK TO THE FLAG LIST
C     IBOUN  = THE LINK FROM THE ENTITY TO THE FLAGS
C     MDIM   = THE DIMENSIONING PARAMETER FOR THE LIST
C
C***********************************************************************
C
      DIMENSION IFLAG (MDIM), IFLIST (MDIM), NEPS (MDIM), LIST (2, MDIM)
      DIMENSION LINKF (2, MDIM), IBOUN (MDIM), LINKE (2, MDIM)
C
      LOGICAL ADDLNK
C
      CALL LTSORT (MDIM, LINKE, JJ, L, ADDLNK)
      IF (L .LE. 0) THEN
         CALL MESAGE ('BOUNDARY CONDITION LINK ATTEMPTED')
         CALL MESAGE ('TO A NONEXISTANT ENTITY')
         WRITE ( * , 10000)IFLAG1, JJ
         CALL MESAGE ('CHECK POINBC,  LINEBC,  OR SIDEBC (S)')
         CALL MESAGE (' ')
      ELSEIF (IBOUN (L) .LE. 0) THEN
         IBOUN (L) = IFLAG1
      ELSE
         K = IBOUN (L)
  100    CONTINUE
         CALL LTSORT (MDIM, LINKF, K, IPNTR, ADDLNK)
         IFLAG2 = IFLAG (IPNTR)
         IF (IFLAG2 .EQ. IFLAG1)GOTO 120
         M1 = IFLIST (IPNTR)
         M2 = M1 + NEPS (IPNTR) - 1
         DO 110 M = M1, M2
            IF (LIST (1, M) .EQ. JJ) THEN
               IF (LIST (2, M) .GT. 0) THEN
                  K = LIST (2, M)
                  GOTO 100
               ELSE
                  LIST (2, M) = IFLAG1
                  GOTO 120
               ENDIF
            ENDIF
  110    CONTINUE
         CALL MESAGE ('PROBLEM LINKING BOUNDARY FLAG TABLES')
      ENDIF
  120 CONTINUE
      RETURN
C
10000 FORMAT (' FLAG: ', I5, '  ENTITY ATTEMPTED:', I5)
C
      END
