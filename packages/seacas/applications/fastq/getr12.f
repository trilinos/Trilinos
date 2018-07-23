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

C $Id: getr12.f,v 1.1 1990/11/30 11:08:37 gdsjaar Exp $
C $Log: getr12.f,v $
C Revision 1.1  1990/11/30 11:08:37  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]GETR12.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GETR12 (MCOM, ICOM, JCOM, CIN, RIN, KIN, R1, R2,
     &   IFOUND)
C***********************************************************************
C
C  SUBROUTINE GETR12 = GETS TWO REAL NUMBERS
C
C***********************************************************************
C
      DIMENSION RIN(MCOM), KIN(MCOM)
      CHARACTER*72 CIN(MCOM)
C
      IF ( (ICOM .GT. JCOM) .OR. (CIN(ICOM) (1:1) .EQ. ' ') ) THEN
         IFOUND = 0
         ICOM = ICOM+1
         R1 = 0.
         R2 = 0.
      ELSE
         R1 = RIN(ICOM)
         ICOM = ICOM+1
         IF ( (ICOM .LE. JCOM) .AND. (KIN(ICOM) .GT. 0) ) THEN
            R2 = RIN(ICOM)
            ICOM = ICOM+1
            IFOUND = 2
         ELSE
            R2 = 0.
            IFOUND = 1
         ENDIF
      ENDIF
      RETURN
C
      END
