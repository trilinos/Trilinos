C Copyright (C) 2009-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C 
C     * Neither the name of NTESS nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C 

C $Id: pltini.f,v 1.2 2000/10/25 13:32:35 gdsjaar Exp $ 
C $Log: pltini.f,v $
C Revision 1.2  2000/10/25 13:32:35  gdsjaar
C Modified intrinsic functions to use generic versions to avoid warnings on SGI 64-bit compiles
C
C Revision 1.1  1993/07/16 16:48:24  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      SUBROUTINE PLTINI(MIN,MAX,START,REND,INTER,EXP,NMIN)
      REAL MIN,MAX,INTER
      INTEGER EXP

      DELTA = MAX - MIN
      IF (DELTA.LT.0.) THEN
         CALL PLTFLU
         CALL SIORPT('PLTINI',
     *               'Maximum value must be greater than minimum value.'
     *               ,2)
         RETURN

      END IF

      CALL PLTINO(MIN,MAX,START,REND,INTER,EXP,NMIN)
      TENEXP = 10.**EXP
      IEXP = NINT(LOG10(ABS(INTER))) - 2
      SMALL = 10.**IEXP
      IF (ABS(START-MIN/TENEXP).GT.SMALL) THEN
         START = START + INTER
      END IF

      IF (ABS(REND-MAX/TENEXP).GT.SMALL) THEN
         REND = REND - INTER
      END IF

      RETURN

      END
