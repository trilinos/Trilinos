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

C $Id: putcrs.f,v 1.2 2000/11/13 15:39:05 gdsjaar Exp $
C $Log: putcrs.f,v $
C Revision 1.2  2000/11/13 15:39:05  gdsjaar
C Cleaned up unused variables and labels.
C
C Removed some real to int conversion warnings.
C
C Revision 1.1.1.1  1990/11/30 11:13:57  gdsjaar
C FASTQ Version 2.0X
C
c Revision 1.1  90/11/30  11:13:56  gdsjaar
c Initial revision
c
C
CC* FILE: [.MAIN]PUTCRS.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE PUTCRS (X, Y, OLDCUR)
C***********************************************************************
C
C     SUBROUTINE PUTCRS = PLACES THE CROSSHAIRS AT THE CURRENT LOCATION
C
C***********************************************************************
C
C
      DIMENSION IDUM(2)
C
      LOGICAL OLDCUR
C
      CHARACTER DUMMY*16
C
C  SELECT DECIMAL MODE
C
      DUMMY = CHAR(27)
      DUMMY(2:4) = 'OR1'
      WRITE(*,*)DUMMY
C
C  PLACE THE CROSSHAIRS AT THE RIGHT LOCATION
C
      CALL MP2PT(1, X, Y, X1, Y1, IDUM)
      IX = INT(X1*4151.)
      IY = INT(Y1*4151.)
      DUMMY(1:1) = CHAR(27)
      DUMMY(2:2) = 'P'
      WRITE(DUMMY(3:8), '(I6)')IX
      DUMMY(9:9) = ','
      WRITE(DUMMY(10:15), '(I6)')IY
      DUMMY(16:16) = ','
      WRITE(*,*)DUMMY
C
C  UNSELECT DECIMAL MODE
C
      DUMMY = CHAR(27)
      DUMMY(2:4) = 'OR0'
      WRITE(*,*)DUMMY
C
      IF(.NOT.OLDCUR)THEN
C
C  ACTIVATE THE CROSSHAIRS
C
         DUMMY = CHAR(27)
         DUMMY(2:3) = 'G1'
         WRITE(*,*)DUMMY
         OLDCUR = .TRUE.
      ENDIF
C
      WRITE(*, '(A)')' '//CHAR(27)//'[2J'
      RETURN
C
      END
