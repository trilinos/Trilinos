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

C $Id: dlpara.f,v 1.1 1990/11/30 11:06:13 gdsjaar Exp $
C $Log: dlpara.f,v $
C Revision 1.1  1990/11/30 11:06:13  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]DLPARA.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE DLPARA (X1, Y1, X2, Y2, XM, B, BAD)
C***********************************************************************
C
C  SUBROUTINE DLPARA = DETERMINES LINE PARAMETERS FROM TWO POINTS
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     INREGN = INPUTS REGION CONNECTIVITIES
C
C***********************************************************************
C
C  VARIABLES USED:
C     X1    = X VALUE OF POINT 1
C     X2    = X VALUE OF POINT 2
C     Y1    = Y VALUE OF POINT 1
C     Y2    = Y VALUE OF POINT 2
C     XM    = THE SLOPE OF A STRIGHT LINE BETWEEN POINT 1 AND 2
C     B     = THE Y INTERCEPT OF THE STRAIGHT LINE BETWEEN POINT 1 AND 2
C
C***********************************************************************
C
      LOGICAL BAD
C
      IF (ABS (X2 - X1) .LT. 0.000001) THEN
         BAD = .TRUE.
         B = X1
      ELSE
         BAD = .FALSE.
         XM =  (Y2 - Y1) / (X2 - X1)
         B = Y1- (X1 * XM)
      ENDIF
      RETURN
      END
