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

C $Id: donrm3.f,v 1.1 1991/02/21 15:42:57 gdsjaar Exp $
C $Log: donrm3.f,v $
C Revision 1.1  1991/02/21 15:42:57  gdsjaar
C Initial revision
C
      SUBROUTINE DONRM3 (COORD, LTNESS, MAP, DIRCOS, TEMP,
     *    NSEG, NUMNIQ, NUMNP)
      DIMENSION COORD(NUMNP, *), LTNESS(4,*), MAP(*), DIRCOS(5,*),
     *    TEMP(3,*)
C
      DO 10 ISEG = 1, NSEG
          XI = COORD(LTNESS(1,ISEG),1 )
          YI = COORD(LTNESS(1,ISEG),2 )
          ZI = COORD(LTNESS(1,ISEG),3 )
C
          XJ = COORD(LTNESS(2,ISEG),1 )
          YJ = COORD(LTNESS(2,ISEG),2 )
          ZJ = COORD(LTNESS(2,ISEG),3 )
C
          XK = COORD(LTNESS(3,ISEG),1 )
          YK = COORD(LTNESS(3,ISEG),2 )
          ZK = COORD(LTNESS(3,ISEG),3 )
C
          XL = COORD(LTNESS(4,ISEG),1 )
          YL = COORD(LTNESS(4,ISEG),2 )
          ZL = COORD(LTNESS(4,ISEG),3 )
C
          AI =  (YK - YI) * (ZL - ZJ) - (ZK - ZI) * (YL - YJ)
          BJ =  (ZK - ZI) * (XL - XJ) - (XK - XI) * (ZL - ZJ)
          CK =  (XK - XI) * (YL - YJ) - (YK - YI) * (XL - XJ)
          RMAG = SQRT ( AI**2 + BJ**2 + CK**2)
C
          TEMP(1,ISEG) = AI / RMAG
          TEMP(2,ISEG) = BJ / RMAG
          TEMP(3,ISEG) = CK / RMAG
C
   10 CONTINUE
C
      DO 20 I=1,NUMNIQ
          DIRCOS(1,I) = 0.0
          DIRCOS(2,I) = 0.0
          DIRCOS(3,I) = 0.0
   20 CONTINUE
C
      DO 40 ISEG = 1, NSEG
          DO 30 J = 1, 4
              MISEG = MAP( 4 * (ISEG-1) + J )
              DIRCOS(1,MISEG) = DIRCOS(1,MISEG) + TEMP(1,ISEG)
              DIRCOS(2,MISEG) = DIRCOS(2,MISEG) + TEMP(2,ISEG)
              DIRCOS(3,MISEG) = DIRCOS(3,MISEG) + TEMP(3,ISEG)
   30     CONTINUE
   40 CONTINUE
C
C ... NORMALIZE ALL DIRECTION COSINES
C
      DO 50 I = 1, NUMNIQ
          A = DIRCOS(1,I)
          B = DIRCOS(2,I)
          C = DIRCOS(3,I)
          R = SQRT(A**2 + B**2 + C**2)
          DIRCOS(1,I) = A/R
          DIRCOS(2,I) = B/R
          DIRCOS(3,I) = C/R
   50 CONTINUE
      RETURN
      END
