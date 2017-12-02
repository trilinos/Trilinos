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

C $Id: indx.f,v 1.1 1990/11/30 11:09:32 gdsjaar Exp $
C $Log: indx.f,v $
C Revision 1.1  1990/11/30 11:09:32  gdsjaar
C Initial revision
C
C
CC* FILE: [.RENUM]INDX.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      FUNCTION INDX (N, L, IVAL)
C************************************************************************
C
C  FUNCTION INDX = FINDS THE INDEX IN L OF IVAL
C
C************************************************************************
C
C  NOTE:
C     L MUST BE IN INCREASING ORDER
C     IF IVAL IS NOT IN L,  INDEX=0 IS RETURNED
C
C***********************************************************************
C
      DIMENSION L (N)
C
C  BISECTION SEARCH
C
      IF (N .LT. 1) THEN
         INDX=0
         RETURN
      ENDIF
      ILO=1
      IHI=N
  100 CONTINUE
      IMID= (ILO + IHI) / 2
C
C  CONVERGENCE
C
      IF (IMID .EQ. ILO) THEN
         IF (IVAL .EQ. L (IMID)) THEN
            INDX=IMID
            RETURN
         ELSEIF (IVAL .NE. L (IHI)) THEN
            INDX=0
            RETURN
         ENDIF
         INDX=IHI
         RETURN
      ENDIF
C
      IF (IVAL .LT. L (IMID)) THEN
         IHI=IMID
      ELSEIF (IVAL .EQ. L (IMID)) THEN
         INDX=IMID
         RETURN
      ELSE
         ILO=IMID
      ENDIF
      GOTO 100
C
      END
