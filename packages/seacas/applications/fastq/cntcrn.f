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

C $Id: cntcrn.f,v 1.2 1999/06/21 22:43:40 gdsjaar Exp $
C $Log: cntcrn.f,v $
C Revision 1.2  1999/06/21 22:43:40  gdsjaar
C Fixed more uninitialized variables; one was causing core dump on g77
C compiled executable.
C
C VERSN was not consistently defined -- now 10 characters everywhere
C
C Updated so full version string output
C
C Added capability to debug memory using unit specified in EXT99
C variable. Similar to STRTUP in SUPLIB
C
C Cleaned up some other code
C
C Upped version
C
C Revision 1.1.1.1  1990/11/30 11:05:07  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:05:06  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]CNTCRN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CNTCRN (MXND, MXCORN, MLN, LNODES, LCORN, NCORN,
     &   NLOOP, N1, ERR)
C***********************************************************************
C
C  SUBROUTINE CNTCRN = COUNTS THE CURRENT DESIGNATED CORNER LENGTHS
C
C***********************************************************************
C
      DIMENSION LNODES (MLN, MXND), LCORN (MXCORN)
C
      LOGICAL ERR
C
      ERR = .FALSE.
C
C  COUNT THE CURRENT CORNERS STARTING AT THE I COUNTER
C
      NODE = N1
      NOLD = N1
      KOUNT = 0
      NCORN = 0
      KKC = 0
      KOUNTC = 0
      LASTC = 0
  100 CONTINUE
      KOUNT = KOUNT + 1
      IF (KOUNT .GT. NLOOP) THEN
         CALL MESAGE ('PROBLEM IN CNTCRN WITH UNCLOSED LOOP')
         ERR = .TRUE.
         GOTO 110
      ENDIF
C
C  A NEW CORNER NODE HAS BEEN FOUND
C
      IF (LNODES (1, NODE) .EQ. 1) THEN
C
C  ADD UP THE NUMBER OF NODES FROM THE LAST "NON-SIDE"
C
         NCORN = NCORN + 1
         IF (NCORN .LE. MXCORN) LCORN(NCORN) = NODE
         IF (NCORN .GT. 1) THEN
            LNODES (7, LASTC) = KOUNTC + 1
         ELSE
            KKC = KOUNTC + 1
         ENDIF
         LASTC = NODE
         KOUNTC = 0
C
C  THIS IS A SIDE - JUST CONTINUE
C
      ELSE
         KOUNTC = KOUNTC + 1
      ENDIF
C
C  CHECK FOR COMPLETION OF THE LOOP
C
      NODE = LNODES (3, NODE)
      IF (NODE .NE. NOLD) GOTO 100
C
C  GET THE FIRST CORNER'S DISTANCE FROM PREVIOUS CORNER CORRECT
C
      LNODES (7, LASTC) = KKC + KOUNTC
C
  110 CONTINUE
C
      RETURN
C
      END
