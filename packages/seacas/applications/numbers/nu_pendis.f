C    Copyright(C) 1988-2017 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
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

C $Id: pendis.f,v 1.1 1991/02/21 15:44:47 gdsjaar Exp $
C $Log: pendis.f,v $
C Revision 1.1  1991/02/21 15:44:47  gdsjaar
C Initial revision
C
      SUBROUTINE PENDIS (SCORD, FCORD, DIST, NDIM, NNOD)
      DIMENSION SCORD(*), FCORD(NDIM, NNOD)

      IF (NDIM .EQ. 3) THEN

C -- DETERMINE PLANE EQUATION

         XI = FCORD(1, 1)
         YI = FCORD(2, 1)
         ZI = FCORD(3, 1)
C
         XJ = FCORD(1, 2)
         YJ = FCORD(2, 2)
         ZJ = FCORD(3, 2)
C
         XK = FCORD(1, 3)
         YK = FCORD(2, 3)
         ZK = FCORD(3, 3)
C
         XL = FCORD(1, 4)
         YL = FCORD(2, 4)
         ZL = FCORD(3, 4)
C
         A =  (YK - YI) * (ZL - ZJ) - (ZK - ZI) * (YL - YJ)
         B =  (ZK - ZI) * (XL - XJ) - (XK - XI) * (ZL - ZJ)
         C =  (XK - XI) * (YL - YJ) - (YK - YI) * (XL - XJ)
         RMAG = SQRT (A**2 + B**2 + C**2)
C
         A = A / RMAG
         B = B / RMAG
         C = C / RMAG
         D = A * FCORD(1,1) + B * FCORD(2,1) + C * FCORD(3,1)

         DIST = ABS(A * SCORD(1) + B * SCORD(2) + C * SCORD(3) - D) /
     *      SQRT(A**2 + B**2 + C**2)

      ELSE IF (NDIM .EQ. 2) THEN
         A  = FCORD(1,2) - FCORD(1,1)
         B  = FCORD(2,2) - FCORD(2,1)

         X1 = FCORD(1,1)
         Y1 = FCORD(2,1)
C
         X0 = SCORD(1)
         Y0 = SCORD(2)
         T  = -1. * (A * (X1 - X0) + B * (Y1 - Y0)) / (A**2 + B**2)
C
         X = X1 + A * T
         Y = Y1 + B * T
C
         DIST = SQRT((X - X0)**2 + (Y - Y0)**2)
      END IF
      RETURN
      END
