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

C $Id: mppola.f,v 1.1 1993/07/16 16:47:23 gdsjaar Exp $
C $Log: mppola.f,v $
C Revision 1.1  1993/07/16 16:47:23  gdsjaar
C Changed plt to library rather than single source file.
C
C=======================================================================
      SUBROUTINE MPPOLA(DIST,AZIM,INC,TWIST)
      COMMON /MAP/MODEL(4,4),VIEW(4,4),PROJ(4,4),CPNEAR,CPFAR,VWPORT(4),
     *       MVP(4,4),VP(4,4),CPLINE(2,2,10),CPPLAN(2,3,10),PEYE(3),
     *       PLOOK(3),ETWIST,NCPLIN,NCPLAN,TMAT1(4,4),TMAT2(4,4),
     *       TMAT3(4,4),TVEC1(4),TVEC2(4),TVEC3(4),TVEC4(4),TARR1(32),
     *       TARR2(32),TARR3(32),TARR4(32),TARR5(32),TARR6(32),
     *       TARR7(32),TARR8(32)
      REAL MODEL,MVP
      REAL INC
      PARAMETER (DPR=57.2958)

      ANG = -AZIM/DPR
      CALL LDROTA('y',COS(ANG),SIN(ANG),TMAT1)
      ANG = -INC/DPR
      CALL LDROTA('x',COS(ANG),SIN(ANG),TMAT2)
      CALL MXMULT(4,TMAT1,TMAT2,TMAT3)
      ANG = -TWIST/DPR
      CALL LDROTA('z',COS(ANG),SIN(ANG),TMAT1)
      CALL MXMULT(4,TMAT3,TMAT1,TMAT2)
      CALL LDTRAN(0.,0.,-DIST,TMAT1)
      CALL MXMULT(4,TMAT2,TMAT1,VIEW)
      CALL MXMULT(4,VIEW,PROJ,VP)
      CALL MXMULT(4,MODEL,VP,MVP)
      RINC = INC/DPR
      RAZIM = AZIM/DPR
      PEYE(1) = DIST*COS(RINC)*SIN(RAZIM)
      PEYE(2) = -DIST*SIN(RINC)
      PEYE(3) = DIST*COS(RINC)*COS(RAZIM)
      PLOOK(1) = 0.
      PLOOK(2) = 0.
      PLOOK(3) = 0.
      ETWIST = TWIST
      RETURN

      END
