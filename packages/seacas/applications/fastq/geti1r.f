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

C $Id: geti1r.f,v 1.1 1990/11/30 11:08:11 gdsjaar Exp $
C $Log: geti1r.f,v $
C Revision 1.1  1990/11/30 11:08:11  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]GETI1R.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GETI1R (MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN, I1, R1,
     &   IFOUND)
C***********************************************************************
C
C  SUBROUTINE GETI1R = GETS AN INTEGER AND A REAL INPUT NUMBER
C
C***********************************************************************
C
      DIMENSION IIN(MCOM), RIN(MCOM), KIN(MCOM)
      CHARACTER*72 CIN(MCOM)
C
      IF( (ICOM .GT. JCOM) .OR. (CIN(ICOM) (1:1) .EQ. ' ') ) THEN
         ICOM = ICOM+1
         IFOUND = 0
      ELSE
         I1 = IIN(ICOM)
         ICOM = ICOM+1
         IF ( (ICOM .LE. JCOM) .AND. (KIN(ICOM) .GT. 0) ) THEN
            R1 = RIN(ICOM)
            ICOM = ICOM + 1
            IFOUND = 2
         ELSE
            R1 = 0.
            IFOUND = 1
         ENDIF
      ENDIF
      RETURN
C
      END
