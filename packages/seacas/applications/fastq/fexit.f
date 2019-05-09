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

C $Id: fexit.f,v 1.3 1999/01/27 15:17:47 gdsjaar Exp $
C $Log: fexit.f,v $
C Revision 1.3  1999/01/27 15:17:47  gdsjaar
C Added typical summary of mesh data on output.
C
C Better filename handling
C
C Cleaned up some character string handling
C
C Revision 1.2  1993/07/21 18:11:37  gdsjaar
C Removed message after stop
C
c Revision 1.1.1.1  1990/11/30  11:07:22  gdsjaar
c FASTQ Version 2.0X
c
c Revision 1.1  90/11/30  11:07:21  gdsjaar
c Initial revision
c
C
CC* FILE: [.MAIN]FEXIT.FOR
CC* MODIFIED BY: TED BLACKER
CC* MODIFICATION DATE: 7/6/90
CC* MODIFICATION: COMPLETED HEADER INFORMATION
C
      SUBROUTINE FEXIT (WROTE, MCOM, ICOM, JCOM, CIN, IIN, RIN, KIN,
     &   TIME1, BATCH, VERSN)
C***********************************************************************
C
C  SUBROUTINE FEXIT = GRACEFULL FASTQ EXIT
C
C***********************************************************************
C
      CHARACTER*72 CIN (MCOM), VERSN*10, DATE*8, TIME*8
C
      LOGICAL IANS, WROTE, BATCH
C
      DIMENSION KIN (MCOM), IIN (MCOM), RIN (MCOM)
C
      IF (.NOT.WROTE)THEN
         CALL MESAGE (' ')
         CALL MESAGE ('***********************************************')
         CALL MESAGE ('*  WARNING: NO OUTPUT FILE HAS BEEN WRITTEN   *')
         CALL MESAGE ('***********************************************')
         CALL INTRUP ('EXIT ANYWAY', IANS, MCOM, ICOM, JCOM, CIN, IIN,
     &      RIN, KIN)
         IF (.NOT.IANS)RETURN
      ENDIF
      CALL MESAGE (' ')
      CALL EXCPUS (TIME2)
      IF (BATCH)THEN
         CALL MESAGE ('FASTQ COMPLETED SUCCESSFULLY')
         CALL EXDATE (DATE)
         CALL EXTIME (TIME)
         WRITE (*, *)'             DATE: ', DATE
         WRITE (*, *)'             TIME: ', TIME
         WRITE (*, *)'          VERSION: ', VERSN
         WRITE (*, ' (A, I5)')'  CPU SECONDS USED: ', IFIX (TIME2-TIME1)
         CALL MESAGE (' ')
      ELSE
         WRITE (*, ' (A, I5)')' CPU SECONDS USED: ', IFIX (TIME2-TIME1)
         CALL MESAGE ('*--------------------------*')
         CALL MESAGE (' ')
         CALL VDESCP (10003, 0, 0)
         CALL PLTEND
      ENDIF
      stop
      END
