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

C $Id: dstort.f,v 1.2 2000/11/13 15:39:04 gdsjaar Exp $
C $Log: dstort.f,v $
C Revision 1.2  2000/11/13 15:39:04  gdsjaar
C Cleaned up unused variables and labels.
C
C Removed some real to int conversion warnings.
C
C Revision 1.1.1.1  1990/11/30 11:06:27  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:06:25  gdsjaar
c Initial revision
c 
C
CC* FILE: [.PAVING]DSTORT.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE DSTORT (X1, X2, X3, X4, Y1, Y2, Y3, Y4, VALUE)
C***********************************************************************
C
C  SUBROUTINE DSTORT = CALCULATES A DISTORTION METRIC FOR AN ELEMENT
C                    USING THE IDEAS IN THE PAPER BY ODDY, 1988.
C
C***********************************************************************
C
C  SETUP THE JACOBIAN MATRIX
C
      XJ11 = (X1 * .125) + (X2 * .375) - (X3 * .375) - (X4 * .125)
      XJ12 = (Y1 * .125) + (Y2 * .375) - (Y3 * .375) - (Y4 * .125)
      XJ21 = - (X1 * .375) + (X2 * .375) + (X3 * .125) - (X4 * .125)
      XJ22 = - (Y1 * .375) + (Y2 * .375) + (Y3 * .125) - (Y4 * .125)
C
C  NORMALIZE THE JACOBIAN WITH RESPECT TO THE ELEMENT SIZE
C
      DETERM = (XJ11 * XJ22) - (XJ12 * XJ21)
      IF (DETERM .LE. 0.) THEN
         VALUE = 1.0E10
         RETURN
      ENDIF
      FACTOR = 1. / SQRT (DETERM)
      XJ11 = XJ11 * FACTOR
      XJ12 = XJ12 * FACTOR
      XJ21 = XJ21 * FACTOR
      XJ22 = XJ22 * FACTOR
C
C  NOW USE THE SECOND INVARIANT OF GREEN'S STRAIN
C
      C11 = XJ11*XJ11 + XJ21*XJ21
      C12 = XJ11*XJ12 + XJ21*XJ22
      C22 = XJ12*XJ12 + XJ22*XJ22
C
      VALUE = C11**2 + 2.*(C12**2) + C22**2 -
     &   (.5 * (C11+C22)**2 )
      VALUE = AMAX1 (VALUE, 0.)
C
      RETURN
C
      END
