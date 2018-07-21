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

C $Id: donrm2.f,v 1.1 1991/02/21 15:42:55 gdsjaar Exp $
C $Log: donrm2.f,v $
C Revision 1.1  1991/02/21 15:42:55  gdsjaar
C Initial revision
C
      SUBROUTINE DONRM2 (COORD, LTNESS, MAP, DIRCOS, TEMP,
     *    NSEG, NUMNIQ, NUMNP)
      DIMENSION COORD(NUMNP, *), LTNESS(2,*), MAP(*), DIRCOS(4,*),
     *    TEMP(2,*)
C
      DO 10 ISEG = 1, NSEG
          XI = COORD( LTNESS(1,ISEG),1 )
          YI = COORD( LTNESS(1,ISEG),2 )
C
          XJ = COORD( LTNESS(2,ISEG),1 )
          YJ = COORD( LTNESS(2,ISEG),2 )
C
          DX = XI - XJ
          DY = YI - YJ
          RMAG = SQRT ( DX**2 + DY**2)
C
          TEMP(1,ISEG) = -DY / RMAG
          TEMP(2,ISEG) =  DX / RMAG
C
   10 CONTINUE
C
      DO 20 I=1,NUMNIQ
          DIRCOS(1,I) = 0.0
          DIRCOS(2,I) = 0.0
   20 CONTINUE
C
      DO 40 ISEG = 1, NSEG
          DO 30 J = 1, 2
              MISEG = MAP( 2 * (ISEG-1) + J )
              DIRCOS(1,MISEG) = DIRCOS(1,MISEG) + TEMP(1,ISEG)
              DIRCOS(2,MISEG) = DIRCOS(2,MISEG) + TEMP(2,ISEG)
   30     CONTINUE
   40 CONTINUE
C
C ... NORMALIZE ALL DIRECTION COSINES
C
      DO 50 I = 1, NUMNIQ
          A = DIRCOS(1,I)
          B = DIRCOS(2,I)
          R = SQRT(A**2 + B**2)
          DIRCOS(1,I) = A/R
          DIRCOS(2,I) = B/R
   50 CONTINUE
      RETURN
      END
