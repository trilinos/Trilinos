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

C $Id: gtxdum.f,v 1.1 1990/11/30 11:09:08 gdsjaar Exp $
C $Log: gtxdum.f,v $
C Revision 1.1  1990/11/30 11:09:08  gdsjaar
C Initial revision
C
C
CC* FILE: [.MAIN]GTXDUM.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE GTXDUM (X, DUMMY, LEN)
C***********************************************************************
C
C  SUBROUTINE GTXDUM = GETS A REAL INTO A DUMMY CHARACTER STRING
C
C***********************************************************************
C
      CHARACTER*72 DUMMY
C
      DUMMY = ' '
      IF(X .LE. -10000.) THEN
         WRITE (DUMMY(1:10), 10000) X
         LEN = 10
      ELSEIF(X .LE. -1000.) THEN
         WRITE (DUMMY(1:6), 10010) X
         LEN = 6
      ELSEIF(X .LE. -100.) THEN
         WRITE (DUMMY(1:6), 10020) X
         LEN = 6
      ELSEIF(X .LE. -10.) THEN
         WRITE (DUMMY(1:6), 10030) X
         LEN = 6
      ELSEIF(X .LE. -1.) THEN
         WRITE (DUMMY(1:6), 10040) X
         LEN = 6
      ELSEIF(X .LT. 0.) THEN
         WRITE (DUMMY(1:6), 10050) X
         LEN = 6
      ELSEIF(X .GE. 10000.) THEN
         WRITE (DUMMY(1:10), 10060) X
         LEN = 9
      ELSEIF(X .GE. 1000.) THEN
         WRITE (DUMMY(1:5), 10070) X
         LEN = 5
      ELSEIF(X .GE. 100.) THEN
         WRITE (DUMMY(1:5), 10080) X
         LEN = 5
      ELSEIF(X .GE. 10.) THEN
         WRITE (DUMMY(1:5), 10090) X
         LEN = 5
      ELSEIF(X .GE. 1.) THEN
         WRITE (DUMMY(1:5), 10100) X
         LEN = 5
      ELSE
         WRITE (DUMMY(1:5), 10110) X
         LEN = 5
      ENDIF
      RETURN
C
10000 FORMAT (1PE10.3)
10010 FORMAT (F6.0)
10020 FORMAT (F6.1)
10030 FORMAT (F6.2)
10040 FORMAT (F6.3)
10050 FORMAT (F6.4)
10060 FORMAT (1PE9.3)
10070 FORMAT (F5.0)
10080 FORMAT (F5.1)
10090 FORMAT (F5.2)
10100 FORMAT (F5.3)
10110 FORMAT (F5.4)
C
      END
