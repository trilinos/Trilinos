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

C $Id: keep34.f,v 1.1 1990/11/30 11:10:43 gdsjaar Exp $
C $Log: keep34.f,v $
C Revision 1.1  1990/11/30 11:10:43  gdsjaar
C Initial revision
C
C
CC* FILE: [.PAVING]KEEP34.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE KEEP34 (ITEST, LTEST, NBEGIN, NEND, ICHNG)
C***********************************************************************
C
C  SUBROTINE KEEP34 = GETS AN ACCEPTABLE SIDE FOR FILLING TO KEEP A
C                     TRIANGLE VALID OR CHANGING TO A RECTANGLE
C
C***********************************************************************
C
      DIMENSION ITEST (3), LTEST (3)
C
C  MAKE SURE THAT THE NBEGIN STARTS AT ONE OF THE CORNERS
C
      IF (NBEGIN .EQ. ITEST(1)) THEN
         NEND = ITEST(2)
      ELSEIF (NBEGIN .EQ. ITEST(2)) THEN
         NEND = ITEST(3)
      ELSEIF (NBEGIN .EQ. ITEST(3)) THEN
         NEND = ITEST(1)
      ELSE
         NBEGIN = ITEST(1)
         NEND = ITEST(3)
      ENDIF
C
C  FIND THE CORRECT ROW (THIS ALREADY ASSUMES THAT THE
C  SUM OF THE SMALLER TWO IS EQUAL TO THE LARGEST ONE)
C
      MMAX = MAX0 (LTEST(1), LTEST(2), LTEST(3))
      IF (LTEST(1) .EQ. MMAX) THEN
         IF (NBEGIN .EQ. ITEST(1)) THEN
            NEND = ICHNG
         ENDIF
      ELSEIF (LTEST(2) .EQ. MMAX) THEN
         IF (NBEGIN .EQ. ITEST(2)) THEN
            NEND = ICHNG
         ENDIF
      ELSE
         IF (NBEGIN .EQ. ITEST(3)) THEN
            NEND = ICHNG
         ENDIF
      ENDIF
C
      RETURN
C
      END
