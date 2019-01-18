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

C $Id: inqstr.f,v 1.3 2000/11/13 15:39:04 gdsjaar Exp $
C $Log: inqstr.f,v $
C Revision 1.3  2000/11/13 15:39:04  gdsjaar
C Cleaned up unused variables and labels.
C
C Removed some real to int conversion warnings.
C
C Revision 1.2  1991/03/22 19:38:52  gdsjaar
C Fixed typo 0 was K0
C
c Revision 1.1.1.1  1990/11/30  11:10:02  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:10:00  gdsjaar
c Initial revision
c 
C
CC* FILE: [.MAIN]UTIL.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE INQSTR (PROMPT, IANS)
C***********************************************************************
C
C  SUBROUTINE INQSTR = INPUTS CHARACTER STRINGS
C
C***********************************************************************
C
      CHARACTER* (*) PROMPT, IANS, HOLD*80
C
      IZ = 0
  100 CONTINUE
      CALL GETINP (IZ, IZ, PROMPT, HOLD, IOSTAT)
      IF (IOSTAT .EQ. 0) THEN
         CALL STRCUT (HOLD)
         IANS = HOLD (1:)
         RETURN
      ELSEIF (IOSTAT .LT. 0) THEN
         IANS = ' '
         RETURN
      ELSEIF (IOSTAT .GT. 0) THEN
         WRITE (*, 10010)
         GOTO 100
      ENDIF
C
10010 FORMAT (' BAD CHARACTER STRING  -  TRY AGAIN')
      END
