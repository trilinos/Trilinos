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

C $Log: tpchk.f,v $
C Revision 1.3  2009/03/25 12:36:48  gdsjaar
C Add copyright and license notice to all files.
C Permission to assert copyright has been granted; blot is now open source, BSD
C
C Revision 1.2  2009/01/22 21:34:22  gdsjaar
C There were several inline dbnums common blocks. Replaced with the
C include so they all have the same size with the added variable types.
C
C Added minor support for nodeset and sideset variables.
C
C It can print the count and the names, but that is all
C at this time.
C
C Revision 1.1  1994/04/07 20:16:27  gdsjaar
C Initial checkin of ACCESS/graphics/blotII2
C
c Revision 1.2  1990/12/14  08:58:58  gdsjaar
c Added RCS Id and Log to all files
c
C=======================================================================
      SUBROUTINE TPCHK (PRTOK, OKAY)
C=======================================================================

C   --*** TPCHK *** (TPLOT) Check that TPLOT is a valid program to run
C   --   Written by Amy Gilkey - revised 11/03/87
C   --
C   --TPCHK checks that TPLOT is a valid program.  This is true if
C   --there are history, global, nodal or element variables in the database.
C   --
C   --If the program is valid, a welcome message is printed.
C   --
C   --Parameters:
C   --   PRTOK - IN - if true, print error messages and welcome messages,
C   --      if false, check but do not print anything
C   --   OKAY - OUT - true iff the program is valid
C   --
C   --Common Variables:
C   --   Uses NVARHI, NVARGL, NVARNP, NVAREL, NSTEPS, NSTEPW of /DBNUMS/

      include 'dbnums.blk'

      LOGICAL PRTOK
      LOGICAL OKAY

C   --Check header information from database

      IF ((NVARHI + NVARGL + NVARNP + NVAREL) .LE. 0) THEN
         IF (PRTOK) CALL PRTERR ('CMDERR',
     &      'No variables are defined')
         GOTO 100
      END IF
      IF ((NSTEPS .LE. 0) .OR.
     &   ((NSTEPW .LE. 0) .AND. (NVARHI .LE. 0))) THEN
         IF (PRTOK) CALL PRTERR ('CMDERR',
     &      'No time steps with variables are defined')
         GOTO 100
      END IF

C   --Write greeting message
      IF (PRTOK) WRITE (*, 10000)

      OKAY = .TRUE.
      RETURN

  100 CONTINUE
      OKAY = .FALSE.
      RETURN
10000  FORMAT (' TPLOT - a time history or X-Y plot program')
      END
