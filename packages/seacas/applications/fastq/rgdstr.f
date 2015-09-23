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

C $Id: rgdstr.f,v 1.2 1998/07/14 18:19:58 gdsjaar Exp $
C $Log: rgdstr.f,v $
C Revision 1.2  1998/07/14 18:19:58  gdsjaar
C Removed unused variables, cleaned up a little.
C
C Changed BLUE labels to GREEN to help visibility on black background
C (indirectly requested by a couple users)
C
C Revision 1.1.1.1  1990/11/30 11:14:57  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:14:56  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]RGDSTR.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE RGDSTR (NPNODE, NPELEM, KKK, NNXK, XN, YN, NXK)
C************************************************************************
C
C  SUBROUTINE RGDSTR = CALCULATES A REGION DISTORTION MEASURE
C
C***********************************************************************
C
      DIMENSION XN (NPNODE), YN (NPNODE), NXK (NNXK, NPELEM)
C
C  CALCULATE THE ELEMENT DISTORTION
C
      N1 = NXK (1,1)
      N2 = NXK (2,1)
      N3 = NXK (3,1)
      N4 = NXK (4,1)
      CALL DSTORT (XN (N1), XN (N2), XN (N3), XN (N4),
     &   YN (N1), YN (N2), YN (N3), YN (N4), VALUE)
      VMIN = VALUE
      VMAX = VALUE
      SUM = VALUE
      KMIN = 1
      KMAX = 1
      DO 100 I = 1, KKK
         N1 = NXK (1,I)
         N2 = NXK (2,I)
         N3 = NXK (3,I)
         N4 = NXK (4,I)
         CALL DSTORT (XN (N1), XN (N2), XN (N3), XN (N4),
     &      YN (N1), YN (N2), YN (N3), YN (N4), VALUE)
         IF (VMIN .GT. VALUE) THEN
            VMIN = VALUE
            KMIN = I
         ELSE IF (VMAX .LT. VALUE) THEN
            VMAX = VALUE
            KMAX = I
         ENDIF
         SUM = SUM + VALUE
  100 CONTINUE
C
C  PRINT OUT THE RESULTS
C
      SUM = SUM / FLOAT (KKK)
      WRITE (*, 10000) VMIN, KMIN, VMAX, KMAX, SUM
C
      RETURN
C
10000 FORMAT ('  THE MINIMUM DISTORTION IS: ',G14.7,' IN ELEMENT: ',I10,
     &   /,
     &   '  THE MAXIMUM DISTORTION IS: ',G14.7,' IN ELEMENT: ',I10, /,
     &   '  THE AVERAGE DISTORTION IS: ',G14.7)
C
      END
