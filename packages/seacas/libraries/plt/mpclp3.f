C Copyright (C) 2009-2017 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
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
C     * Neither the name of NTESS nor the names of its
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

C $Id: mpclp3.f,v 1.1 1993/07/16 16:47:08 gdsjaar Exp $ 
C $Log: mpclp3.f,v $
C Revision 1.1  1993/07/16 16:47:08  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      LOGICAL FUNCTION MPCLP3(N,PX,PY,PZ,VX,VY,VZ)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      CHARACTER*6 SUBNAM
      DIMENSION PX(*),PY(*),PZ(*),VX(*),VY(*),VZ(*)
      PARAMETER (SUBNAM='MPCLP3')

      MPCLP3 = .FALSE.
      IF (N.GT.10) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *               'Too many clipping planes specified; max is 10',2)
         RETURN

      END IF

      IF (N.LT.0) THEN
         CALL PLTFLU
         CALL SIORPT(SUBNAM,
     *               'You cannot specify less than zero clipping planes'
     *               ,2)
         RETURN

      END IF

      MPCLP3 = .TRUE.
      NCPLAN = N
      DO 2440 I = 1,N
         CPPLAN(1,1,I) = PX(I)
         CPPLAN(1,2,I) = PY(I)
         CPPLAN(1,3,I) = PZ(I)
         CPPLAN(2,1,I) = VX(I)
         CPPLAN(2,2,I) = VY(I)
         CPPLAN(2,3,I) = VZ(I)
 2440 CONTINUE
      RETURN

      END
