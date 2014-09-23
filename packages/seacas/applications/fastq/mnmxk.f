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

C $Id: mnmxk.f,v 1.1 1990/11/30 11:12:22 gdsjaar Exp $
C $Log: mnmxk.f,v $
C Revision 1.1  1990/11/30 11:12:22  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]MNMXK.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE MNMXK (NPELEM, NPNODE, NNXK, NXK, XN, YN, CENTK, KKK,
     &   XMIN, XMAX, YMIN, YMAX)
C***********************************************************************
C
C  SUBROUTINE MNMXK = FINDS MIN AND MAX DIMENSIONS FOR FLAGGED ELEMENTS
C
C**********************************************************************
C
      DIMENSION NXK (NNXK, NPELEM), CENTK (2, NPELEM)
      DIMENSION XN (NPNODE), YN (NPNODE)
C
C  FIND THE FIRST ELEMENT TO BE PLOTTED
C
      DO 150 I = 1, KKK
         IF (NXK (1, I) .GT. 0) THEN
            JX1 = I
            JX2 = I
            JY1 = I
            JY2 = I
C
C  COMPARE CENTERS TO GET MIN AND MAX ELEMENTS
C
            DO 100 J = I + 1, KKK
               IF (NXK (1, J) .GT. 0) THEN
                  IF (CENTK (1, J) .LT. CENTK (1, JX1))JX1 = J
                  IF (CENTK (1, J) .GT. CENTK (1, JX2))JX2 = J
                  IF (CENTK (2, J) .LT. CENTK (2, JY1))JY1 = J
                  IF (CENTK (2, J) .GT. CENTK (2, JY2))JY2 = J
               ENDIF
  100       CONTINUE
C
C  FIND CORRECT MIN AND MAX FROM NODES OF MIN AND MAX ELEMENTS
C
            XMIN = XN (NXK (1, JX1))
            DO 110 K = 2, NNXK
               IF (NXK (K, JX1) .GT. 0)
     &            XMIN = AMIN1 (XMIN, XN (NXK (K, JX1)))
  110       CONTINUE
            XMAX = XN (NXK (1, JX2))
            DO 120 K = 2, NNXK
               IF (NXK (K, JX2) .GT. 0)
     &            XMAX = AMAX1 (XMAX, XN (NXK (K, JX2)))
  120       CONTINUE
            YMIN = YN (NXK (1, JY1))
            DO 130 K = 2, NNXK
               IF (NXK (K, JY1) .GT. 0)
     &            YMIN = AMIN1 (YMIN, YN (NXK (K, JY1)))
  130       CONTINUE
            YMAX = YN (NXK (1, JY2))
            DO 140 K = 2, NNXK
               IF (NXK (K, JY2) .GT. 0)
     &            YMAX = AMAX1 (YMAX, YN (NXK (K, JY2)))
  140       CONTINUE
            RETURN
         ENDIF
  150 CONTINUE
C
C  RETURN WITH DEFAULT MINS AND MAXS
C
      XMIN = 0.
      XMAX = 1.
      YMIN = 0.
      YMAX = 1.
      RETURN
C
      END
