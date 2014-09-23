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

C $Id: chcond.f,v 1.1 1990/11/30 11:04:25 gdsjaar Exp $
C $Log: chcond.f,v $
C Revision 1.1  1990/11/30 11:04:25  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]CHCOND.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CHCOND (NPER, NSA, SMANG, INDEX, IFIRST, N1, N2, N3,
     &   N4, CHECK)
C*********************************************************************
C
C  SUBROUTINE CHCOND = THIS SUBROUTINE CHECKS IF THE "NSA" ANGLES
C                       SATISFIES THE CONDITIONS.
C
C*********************************************************************
C
C  VARIABLES  IN: NPER .... NUMBER OF POINTS IN THE REGION
C                 NSA ..... NUMBER OF SOUGHT SMALLEST ANGLES
C                 SMANG ... ARRAY OF SMALLEST ANGLES
C               INDEX ... POINTERS INTO THE ANGLE ARRAY
C           OUT: IFIRST... POINTER TO THE FIRST VERTEX
C                 Mi ...... INTERVALS FOR THE PENTAGON REGION
C                 CHECK ... .EQ. TRUE IF IT SATISFIES THE CONDITIONS
C
C  CALL BY: PICKM5.FOR
C
C WRITTEN BY: HORACIO RECALDE                          DATE:FEB 15, 1988
C
C************************************************************************
C
      PARAMETER (NSANG = 10)
      DIMENSION SMANG(NSA + 1), NAUX(NSANG), INDEX(NSA + 1)
      LOGICAL CHECK
C
      NSA2 = NSA/2
C
C--- SORT THE INDEX ARRAY TO FIND THE 'NSA2' SMALLEST ANGLES
C
      CALL SORTIA (NSA, INDEX, NSA2, NAUX)
      IFIRST = NAUX(1)
      N1 = NAUX(2) - NAUX(1)
      N2 = NAUX(3) - NAUX(2)
      N3 = NAUX(4) - NAUX(3)
      N4 = NAUX(5) - NAUX(4)
      N5 = NPER - N1 - N2 - N3 - N4
C
C--- CHECK COMPATIBILITY EQUATIONS
C
      IF ((N1 + N2 + N3 .LT. N4 + N5 + 2) .OR.
     &   (N2 + N3 + N4 .LT. N5 + N1 + 2) .OR.
     &   (N3 + N4 + N5 .LT. N1 + N2 + 2) .OR.
     &   (N4 + N5 + N1 .LT. N2 + N3 + 2) .OR.
     &   (N5 + N1 + N2 .LT. N3 + N4 + 2)) THEN
         CHECK = .FALSE.
      ELSE
         CHECK = .TRUE.
      ENDIF
C
      RETURN
      END
