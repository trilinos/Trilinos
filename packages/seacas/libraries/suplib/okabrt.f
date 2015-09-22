C Copyright(C) 2009 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software.
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
C     * Neither the name of Sandia Corporation nor the names of its
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

C=======================================================================
      LOGICAL FUNCTION OKABRT (ISOK)
C=======================================================================
C$Id: okabrt.f,v 1.3 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: okabrt.f,v $
CRevision 1.3  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.2  1990/11/30 09:51:00  gdsjaar
CModified to work on Unicos
C
c Revision 1.1.1.1  90/08/14  16:15:56  gdsjaar
c Testing
c 
c Revision 1.1  90/08/14  16:15:54  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:41  gdsjaar
c Initial revision
c 

C   --*** OKABRT *** (ETCLIB) Initialize cancel function
C   --   Written by Amy Gilkey - revised 12/21/87
C   --
C   --OKABRT initializes the cancel flag.  It must be called before ISABRT.

C   --Routines Called:
C   --   CPUIFC - (PLTLIB) Check interrupt flag

      LOGICAL ISABRT
      LOGICAL ISOK

      LOGICAL CPUIFC, LDUM

      LOGICAL DOABRT
      SAVE DOABRT

      DATA DOABRT / .FALSE. /

C   --Initialize enable cancel flag
      DOABRT = ISOK

      IF (DOABRT) THEN
C      --Initialize cancel flag
         LDUM = CPUIFC (.TRUE.)
      END IF

      OKABRT = DOABRT

      RETURN

C=======================================================================
      ENTRY ISABRT ()
C=======================================================================
C$Id: okabrt.f,v 1.3 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: okabrt.f,v $
CRevision 1.3  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.2  1990/11/30 09:51:00  gdsjaar
CModified to work on Unicos
C
c Revision 1.1.1.1  90/08/14  16:15:56  gdsjaar
c Testing
c 
c Revision 1.1  90/08/14  16:15:54  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:41  gdsjaar
c Initial revision
c 

C   --*** ISABRT *** (ETCLIB) Check cancel function
C   --   Written by Amy Gilkey - revised 12/17/87
C   --
C   --ISABRT checks the cancel flag.  If it is set, it aborts the current
C   --processing.  In any case, the value of the cancel flag is returned
C   --as the function value.

C   --Routines Called:
C   --   CPUIFC - (PLTLIB) Check interrupt flag

      IF (DOABRT) THEN
C      --Return cancel flag
         ISABRT = CPUIFC (.FALSE.)

         IF (ISABRT) THEN
C         --Send abort message
            WRITE (*, '(1X, A)') '*** Processing aborted ***'
         END IF

      ELSE
         ISABRT = .FALSE.
      END IF

      RETURN
      END
