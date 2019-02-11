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

C $Id: check.f,v 1.1 1990/11/30 11:04:27 gdsjaar Exp $
C $Log: check.f,v $
C Revision 1.1  1990/11/30 11:04:27  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]CHECK.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE CHECK (MIN, MAX, I)
C***********************************************************************
C
C  SUBROUTINE CHECK = CHECKS 2 VALUES FOR BEING OUT OF PRESCRIBED BOUNDS
C
C***********************************************************************
C
C  SUBROUTINE CALLED BY:
C     LIST   = LISTS POINTS,  LINES,  AND REGIONS USED IN MESH DEFINITION
C     ERASE  = DELETES POINTS,  LINES,  AND REGIONS FROM THE MESH
C              DEFINITIONS
C
C***********************************************************************
C
C  VARIABLES USED:
C     MIN   = MINIMUM VALUE TO BE TESTED
C     MAX   = MAXIMUM VALUE TO BE TESTED
C     I     = THE ABSOLUTE MAXIMUM VALUE ALLOWED  (THE MINIMUM IS 1)
C
C************************************************************************
C
      IF (MIN .LT. 1)MIN = 1
      IF (MAX .GT. I)MAX = I
      IF (MAX .LT. MIN)MAX = MIN
      RETURN
      END
