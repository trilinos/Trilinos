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

C $Id: pltcnm.f,v 1.1 1993/07/16 16:47:48 gdsjaar Exp $ 
C $Log: pltcnm.f,v $
C Revision 1.1  1993/07/16 16:47:48  gdsjaar
C Changed plt to library rather than single source file.
C 
C=======================================================================
      LOGICAL FUNCTION PLTCNM(VALUE,CLR)
      CHARACTER*(*) CLR

      PLTCNM = .TRUE.
      IF (VALUE.EQ.0.) THEN
         CLR = 'BLACK'

      ELSE IF (VALUE.EQ.1.) THEN
         CLR = 'RED'

      ELSE IF (VALUE.EQ.2.) THEN
         CLR = 'GREEN'

      ELSE IF (VALUE.EQ.3.) THEN
         CLR = 'YELLOW'

      ELSE IF (VALUE.EQ.4.) THEN
         CLR = 'BLUE'

      ELSE IF (VALUE.EQ.6.) THEN
         CLR = 'CYAN'

      ELSE IF (VALUE.EQ.5.) THEN
         CLR = 'MAGENTA'

      ELSE IF (VALUE.EQ.7.) THEN
         CLR = 'WHITE'

      ELSE IF (VALUE.EQ.8.) THEN
         CLR = 'GRAY'

      ELSE IF (VALUE.EQ.10.) THEN
         CLR = 'DKGRAY'

      ELSE IF (VALUE.EQ.9.) THEN
         CLR = 'LTGRAY'

      ELSE IF (VALUE.EQ.12.) THEN
         CLR = 'LIME'

      ELSE IF (VALUE.EQ.11.) THEN
         CLR = 'PINK'

      ELSE IF (VALUE.EQ.15.) THEN
         CLR = 'ORANGE'

      ELSE IF (VALUE.EQ.14.) THEN
         CLR = 'VIOLET'

      ELSE IF (VALUE.EQ.13.) THEN
         CLR = 'LTBLUE'

      ELSE
         PLTCNM = .FALSE.
         CLR = 'UNKNOWN'
      END IF

      RETURN

      END
