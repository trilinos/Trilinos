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
      SUBROUTINE DBPNAM (OPTION, NVARHI, NVARGL, NVARNP, NVAREL,
     &   NAMEHV, NAMEGV, NAMENV, NAMEEV)
C=======================================================================
C$Id: dbpnam.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: dbpnam.f,v $
CRevision 1.2  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:13:51  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:13:49  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:19  gdsjaar
c Initial revision
c 

C   --*** DBPNAM *** (EXOLIB) Print database variable names
C   --   Written by Amy Gilkey - revised 01/21/88
C   --
C   --DBPNAM displays the database variable names.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options:
C   --      'H' to print history variable names
C   --      'G' to print global variable names
C   --      'N' to print nodal variable names
C   --      'E' to print element variable names
C   --   NVARHI - IN - the number of history variables (if OPTION)
C   --   NVARGL - IN - the number of global variables (if OPTION)
C   --   NVARNP - IN - the number of nodal variables (if OPTION)
C   --   NVAREL - IN - the number of element variables (if OPTION)
C   --   NAMEHV - IN - the history variable names (if OPTION)
C   --   NAMEGV - IN - the global variable names (if OPTION)
C   --   NAMENV - IN - the nodal variable names (if OPTION)
C   --   NAMEEV - IN - the element variable names (if OPTION)

      PARAMETER (NAMSPC = 8+6)

      CHARACTER*(*) OPTION
      INTEGER NVARHI, NVARGL, NVARNP, NVAREL
      CHARACTER*8 NAMEHV(*)
      CHARACTER*8 NAMEGV(*)
      CHARACTER*8 NAMENV(*)
      CHARACTER*8 NAMEEV(*)

      WRITE (*, 10000)

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'H') .GT. 0)) THEN
         WRITE (*, 10010) 'History:', (NAMEHV(I), I=1,NVARHI)
      END IF
      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'G') .GT. 0)) THEN
         WRITE (*, 10010) 'Global: ', (NAMEGV(I), I=1,NVARGL)
      END IF
      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'N') .GT. 0)) THEN
         WRITE (*, 10010) 'Nodal:  ', (NAMENV(I), I=1,NVARNP)
      END IF
      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'E') .GT. 0)) THEN
         WRITE (*, 10010) 'Element:', (NAMEEV(I), I=1,NVAREL)
      END IF

      RETURN

10000  FORMAT (/, 1X, 'Variables Names:')
10010  FORMAT (4X, A, :, 6 (2X, A8), :, /,
     &   (12X, 6 (2X, A8)))
      END
