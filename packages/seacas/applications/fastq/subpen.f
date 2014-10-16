C    Copyright (c) 2014, Sandia Corporation.
C    Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
C    the U.S. Government retains certain rights in this software.
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

C $Id: subpen.f,v 1.1 1990/11/30 11:16:52 gdsjaar Exp $
C $Log: subpen.f,v $
C Revision 1.1  1990/11/30 11:16:52  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]SUBPEN.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE SUBPEN (NPER, NEWPER, X, Y, NID, XSUB, YSUB, NIDSUB,
     &   NUM, M1, M2, IADD, ITRI, XCEN, YCEN)
C***********************************************************************
C
C  SUBROUTINE SUBPEN = PUTS A PENTAGON SUBREGION'S PERIMETER INTO THE
C                      NPERIM ARRAYS
C
C***********************************************************************
C
      DIMENSION X (NPER), Y (NPER), NID (NPER)
      DIMENSION XSUB (NPER), YSUB (NPER), NIDSUB (NPER)
C
C  PUT SIDE ONE AND TWO INTO THE PERIMETER LIST
C
      KOUNT = 0
      DO 100 I = 1, NUM + 1
         KOUNT = KOUNT + 1
         J = I + IADD
         IF (J .GT. NPER)J = J - NPER
         XSUB (KOUNT) = X (J)
         YSUB (KOUNT) = Y (J)
         NIDSUB (KOUNT) = NID (J)
  100 CONTINUE
C
C  PUT SIDE THREE INTO THE LIST
C
      XDIF = XCEN - XSUB (KOUNT)
      YDIF = YCEN - YSUB (KOUNT)
      XINT = XDIF / FLOAT (M1)
      YINT = YDIF / FLOAT (M1)
      DO 110 I = 1, M1 - 1
         KOUNT = KOUNT + 1
         XSUB (KOUNT) = XSUB (KOUNT - 1) + XINT
         YSUB (KOUNT) = YSUB (KOUNT - 1) + YINT
         NIDSUB (KOUNT) =  (ITRI * 100000) + M1 - I + 1
  110 CONTINUE
C
C  ENTER THE CENTER POINT
C
      KOUNT = KOUNT + 1
      XSUB (KOUNT) = XCEN
      YSUB (KOUNT) = YCEN
      NIDSUB (KOUNT) = 100000
C
C  PUT SIDE FOUR INTO THE LIST
C
      ITRI2 = ITRI + 2
      IF (ITRI2 .GT. 3)ITRI2 = ITRI2 - 3
      XDIF = X (IADD + 1) - XCEN
      YDIF = Y (IADD + 1) - YCEN
      XINT = XDIF / FLOAT (M2)
      YINT = YDIF / FLOAT (M2)
      DO 120 I = 1, M2 - 1
         KOUNT = KOUNT + 1
         XSUB (KOUNT) = XSUB (KOUNT - 1) + XINT
         YSUB (KOUNT) = YSUB (KOUNT - 1) + YINT
         NIDSUB (KOUNT) =  (100000 * ITRI2) + I + 1
  120 CONTINUE
C
      NEWPER = KOUNT
C
      RETURN
C
      END
