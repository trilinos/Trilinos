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

C $Id: apalib.f,v 1.1 1990/11/30 11:03:34 gdsjaar Exp $
C $Log: apalib.f,v $
C Revision 1.1  1990/11/30 11:03:34  gdsjaar
C Initial revision
C
C
CC* FILE: [.QMESH]APALIB.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE APALIB (MXND, XN, YN, LXK, NXL, K, NODES, AREA, XCEN,
     &   YCEN, KNUM, KLIB, NLIB, ALIB, XCLIB, YCLIB)
C***********************************************************************
C
C  SUBROUTINE APALIB = LIBRARY OF ELEMENT DATA USED TO AVOID DUPLICATE
C                      COMPUTATIONS
C
C***********************************************************************
C
      DIMENSION NODES (4), KLIB (8), NLIB (4, 8)
      DIMENSION ALIB (8), XCLIB (8), YCLIB (8)
      DIMENSION XN (MXND), YN (MXND), LXK (4, MXND), NXL (2, 3 * MXND)
C
      LOGICAL CCW
C
C  SEARCH LIBRARY
C
      IF (KNUM .GT. 0) THEN
         DO 110 I = 1, KNUM
            IF  (K - KLIB (I) .EQ. 0) THEN
C
C  FETCH FROM LIBRARY
C
               IK = I
               DO 100 J = 1, 4
                  NODES (J) = NLIB (J, IK)
  100          CONTINUE
               AREA = ALIB (IK)
               XCEN = XCLIB (IK)
               YCEN = YCLIB (IK)
               RETURN
            ENDIF
  110    CONTINUE
      ENDIF
C
C  COMPUTE NEW DATA
C
      CCW = .TRUE.
      CALL GNXKA (MXND, XN, YN, K, NODES, AREA, LXK, NXL, CCW)
      N1 = NODES (1)
      N2 = NODES (2)
      N3 = NODES (3)
      N4 = NODES (4)
      XCEN =  (XN (N1) + XN (N2) + XN (N3) + XN (N4)) * 0.25
      YCEN =  (YN (N1) + YN (N2) + YN (N3) + YN (N4)) * 0.25
      IF (KNUM .GE. 8) RETURN
C
C  FILE NEW DATA IN LIBRARY
C
      KNUM = KNUM + 1
      DO 120 I = 1, 4
         NLIB (I, KNUM) = NODES (I)
  120 CONTINUE
      KLIB (KNUM) = K
      ALIB (KNUM) = AREA
      XCLIB (KNUM) = XCEN
      YCLIB (KNUM) = YCEN
      RETURN
C
      END
