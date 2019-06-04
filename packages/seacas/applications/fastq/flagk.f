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

C $Id: flagk.f,v 1.1 1990/11/30 11:07:35 gdsjaar Exp $
C $Log: flagk.f,v $
C Revision 1.1  1990/11/30 11:07:35  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]FLAGK.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE FLAGK (NPELEM, NNXK, NXK, MAPDXG, I1, I2, SETFLG, OLD)
C************************************************************************
C
C  SUBROUTINE FLAGK = FLAGS ELEMENTS FOR PLOTTING OR NOT PLOTTING
C
C***********************************************************************
C
C  VARIABLES USED:
C     NPELEMS = NUMBER OF PROCESSED ELEMENTS
C     NXK     = NODES PER ELEMENT ARRAY (CONNECTIVITY)
C     I1      = BEGINNING ELEMENT TO BE FLAGGED
C     I2      = ENDING ELEMENT TO BE FLAGGED
C     SETFLG  = .TRUE. IF THE ELEMENT IS TO BE FLAGGED FOR PLOTTING
C     OLD     = .TRUE. IF THE OLD ELEMENT NUMBERS ARE TO BE USED
C
C***********************************************************************
C
C  NOTE:
C     THE ELEMENT IS FLAGGED FOR PLOTTING BY FORCING THE FIRST NODE TO
C     BE POSITIVE AND VICE VERSUS
C
C***********************************************************************
C
      DIMENSION NXK(NNXK,NPELEM), MAPDXG(NPELEM)
C
      LOGICAL SETFLG, OLD
C
      IF (OLD) THEN
         IF (SETFLG) THEN
            DO 100 I = I1, I2
               NXK(1,I) = IABS (NXK (1,I))
  100       CONTINUE
         ELSE
            DO 110 I = I1, I2
               NXK(1,I) = -IABS (NXK(1,I))
  110       CONTINUE
         ENDIF
      ELSE
         IF (SETFLG) THEN
            DO 120 I = I1, I2
               J = MAPDXG (I)
               NXK(1,J) = IABS (NXK (1,J))
  120       CONTINUE
         ELSE
            DO 130 I = I1, I2
               J = MAPDXG(I)
               NXK(1,J) = -IABS (NXK (1,J))
  130       CONTINUE
         ENDIF
      ENDIF
      RETURN
C
      END
