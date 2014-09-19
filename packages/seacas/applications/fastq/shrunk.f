C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Governement retains certain rights in this software.
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

C $Id: shrunk.f,v 1.1 1990/11/30 11:15:38 gdsjaar Exp $
C $Log: shrunk.f,v $
C Revision 1.1  1990/11/30 11:15:38  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]SHRUNK.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      LOGICAL FUNCTION SHRUNK (RATIO, NROW)
C***********************************************************************
C
C  FUNCTION SHRUNK = LOGICAL FUNCTION THAT RETURNS TRUE IF THE ELEMENT
C                    SIZE IS DIMINISHING WITH ROW DEPTH
C
C***********************************************************************
C
      DATA TOLER1 /.85/, TOLER2 /.75/, TOLER3 /.6/
C
      IF ((NROW .GE. 3) .AND. (RATIO .LT. TOLER1)) THEN
         SHRUNK = .TRUE.
      ELSEIF ((NROW .GE. 2) .AND. (RATIO .LT. TOLER2)) THEN
         SHRUNK = .TRUE.
      ELSEIF ((NROW .GE. 1) .AND. (RATIO .LT. TOLER3)) THEN
         SHRUNK = .TRUE.
      ELSE
         SHRUNK = .FALSE.
      ENDIF
C
      RETURN
      END
