C Copyright (C) 2009 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C     * Redistributions of source code must retain the above copyright
C       notice, this list of conditions and the following disclaimer.
C 
C     * Redistributions in binary form must reproduce the above
C       copyright notice, this list of conditions and the following
C       disclaimer in the documentation and/or other materials provided
C       with the distribution.
C 
C     * Neither the name of Sandia Corporation nor the names of its
C       contributors may be used to endorse or promote products derived
C       from this software without specific prior written permission.
C 
C THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
C "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
C LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
C A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
C OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
C SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
C LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
C DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
C THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
C (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
C OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
C 

C $Id: mppers.f,v 1.1 1993/07/16 16:47:22 gdsjaar Exp $ 
C $Log: mppers.f,v $
C Revision 1.1  1993/07/16 16:47:22  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      LOGICAL FUNCTION MPPERS(FOVY,ASPECT,NEAR,FAR)
      CHARACTER*6 SUBNAM
      PARAMETER (SUBNAM='MPPERS')
      PARAMETER (DPR=57.2958)
      REAL NEAR
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP

      MPPERS = .FALSE.
      IF (NEAR.LT.0.) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *         'You cannot specify the near clipping plane less than 0.'
     *               ,2)
         RETURN

      END IF

      IF (NEAR.GE.FAR) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *'You cannot specify the near clipping plane >= to the far clipping
     * plane',2)
         RETURN

      END IF

      MPPERS = .TRUE.
      CALL MXZERO(4,PROJ)
      PROJ(2,2) = 1./TAN((FOVY/DPR)/2.)
      PROJ(1,1) = PROJ(2,2)/ASPECT
      PROJ(3,3) = - (FAR+NEAR)/ (FAR-NEAR)
      PROJ(4,3) = - (2.*FAR*NEAR)/ (FAR-NEAR)
      PROJ(3,4) = -1.
      CPNEAR = NEAR
      CPFAR = FAR
      CALL MXMULT(4,VIEW,PROJ,VP)
      CALL MXMULT(4,MODEL,VP,MVP)
      RETURN

      END
