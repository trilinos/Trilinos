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

C $Id: ch3to4.f,v 1.1 1990/11/30 11:04:22 gdsjaar Exp $
C $Log: ch3to4.f,v $
C Revision 1.1  1990/11/30 11:04:22  gdsjaar
C Initial revision
C
C

CC* FILE: [.PAVING]CH3TO4.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CH3TO4 (MXND, MXCORN, MLN, NCORN, LCORN, LNODES, ICOMB,
     &   ANGLE, ITEST, LTEST, QUAL, POSBL4, ICHANG)
C***********************************************************************
C
C  SUBROTINE CH3TO4 = CHECKS THE FEASIBILITY OF A
C                     RECTANGLE FROM A TRIANGLE
C
C***********************************************************************
C
      DIMENSION LNODES (MLN, MXND), ANGLE (MXND), LCORN (MXCORN)
      DIMENSION ICOMB (MXCORN), ITEST (3), LTEST (3)
C
      LOGICAL POSBL4
C
C  ASSUME PERFECT QUALITY
C
C      QUAL = 0.
      POSBL4 = .TRUE.
C
C  FIND THE POSSIBLE RECTANGLE (THIS ALREADY ASSUMES THAT THE
C  SUM OF THE SMALLER TWO IS EQUAL TO THE LARGEST ONE)
C
      MMAX = MAX0 (LTEST(1), LTEST(2), LTEST(3))
      IF (LTEST(1) .EQ. MMAX) THEN
         ICHANG = JUMPLP (MXND, MLN, LNODES, ITEST(1), LTEST(2))
      ELSEIF (LTEST(2) .EQ. MMAX) THEN
         ICHANG = JUMPLP (MXND, MLN, LNODES, ITEST(2), LTEST(3))
      ELSE
         ICHANG = JUMPLP (MXND, MLN, LNODES, ITEST(3), LTEST(1))
      ENDIF
C
C  TEST THE POSSIBLE RECTANGLE FOR GOODNESS
C  ADD UP THE NICKS FOR BAD ANGLES AT THE GIVEN CORNERS
C
C      DO 100 I = 1, NCORN
C         IF (ICOMB (I) .EQ. 1) THEN
C            QUAL = QUAL + (.8 * NICKC (ANGLE (LCORN (I)) ))
C         ELSE
C            QUAL = QUAL + (.8 * NICKS (ANGLE (LCORN (I)) ))
C         ENDIF
C  100 CONTINUE
C
C  ADD UP THE NICKS FOR THE NEW CORNER
C
C      QUAL = QUAL + (.8 * NICKS (ANGLE (ICHANG)) )
C
      RETURN
C
      END
