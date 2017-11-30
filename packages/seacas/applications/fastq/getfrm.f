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

C $Id: getfrm.f,v 1.1 1990/11/30 11:08:05 gdsjaar Exp $
C $Log: getfrm.f,v $
C Revision 1.1  1990/11/30 11:08:05  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]GETFRM.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GETFRM (MXND, LINES, NL, NXL, NODE, N0, N2, NFROM)
C***********************************************************************
C
C  SUBROUTINE GETFRM = GETS THE NODES THAT THE CURRENT NODE CAME FROM
C
C***********************************************************************
C
      DIMENSION NXL(2, 3*MXND), LINES(NL)
C
      NFROM = 0
C
      IF (NL .EQ. 3) THEN
         DO 100 IL = 1, NL
            ILL = LINES (IL)
            IF (NXL (1, ILL) .EQ. NODE) THEN
               NTEST = NXL (2, ILL)
            ELSEIF (NXL (2, ILL) .EQ. NODE) THEN
               NTEST = NXL (1, ILL)
            ELSE
               CALL MESAGE ('** PROBLEMS IN GETFRM **')
               GOTO 110
            ENDIF
            IF ((NTEST .NE. N0) .AND. (NTEST .NE. N2)) THEN
               NFROM = NTEST
               GOTO 110
            ENDIF
  100    CONTINUE
      ENDIF
  110 CONTINUE
C
      RETURN
C
      END
