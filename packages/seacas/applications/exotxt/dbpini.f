C Copyright (c) 2007 Sandia Corporation. Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Governement
C retains certain rights in this software.
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
C 

C$Id: dbpini.f,v 1.4 2007/10/17 18:46:10 gdsjaar Exp $
C$Log: dbpini.f,v $
CRevision 1.4  2007/10/17 18:46:10  gdsjaar
CAdded copyright notice to all files.
C
Cexotxt2 is licensed under the BSD license
C
CRevision 1.3  1996/05/21 16:52:22  caforsy
CAdded read/write for property data.  Cleaned up exodusII error checks
C
CRevision 1.2  1995/11/07 15:01:44  gdsjaar
CInitial checkin of ACCESS/translate/exotxt2
C=======================================================================
      SUBROUTINE DBPINI (OPTION, NDB, TITLE, NDIM, NUMNP, NUMEL, NELBLK,
     &   NUMNPS, LNPSNL, LNPSDF, NUMESS, LESSEL, LESSNL, LESSDF,
     &   NVARGL, NVARNP, NVAREL)
C=======================================================================
C   --*** DBPINI *** (EXOLIB) Display database title and initial variables
C   --   Written by Amy Gilkey - revised 10/15/87
C   -- Revised 8/5/95 to work with EXODUSIIV2 database format
C   --
C   --DBPINI displays the database filename (optional) and title and
C   --the database initial variables.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options:
C   --      'N' to print database name
C   --      'T' to print title
C   --      'I' to print number of nodes, etc.
C   --      'S' to print node set and side set information
C   --      'V' to print number of variables
C   --   NDB    - IN - the database file (if OPTION)
C   --   TITLE  - IN - the database title (if OPTION)
C   --   NDIM   - IN - the number of coordinates per node (if OPTION)
C   --   NUMNP  - IN - the number of nodes (if OPTION)
C   --   NUMEL  - IN - the number of elements (if OPTION)
C   --   NELBLK - IN - the number of element blocks (if OPTION)
C   --   NUMNPS - IN - the number of node sets (if OPTION)
C   --   LNPSNL - IN - the length of the node sets node list (if OPTION)
C   --   LNPSDF - IN - the length of the node set distribution list
C   --   NUMESS - IN - the number of side sets (if OPTION)
C   --   LESSEL - IN - the length of the side sets element list (if OPTION)
C   --   LESSNL - IN - the length of the side sets node list (if OPTION)
C   --   LESSDF - IN - the length of the side set distribution factors list
C   --   NVARGL - IN - the number of global variables (if OPTION)
C   --   NVARNP - IN - the number of nodal variables (if OPTION)
C   --   NVAREL - IN - the number of element variables (if OPTION)

      CHARACTER*80 FILNAM
      CHARACTER*(*) OPTION
      CHARACTER*80 TITLE
      INTEGER NDB
      INTEGER NDIM, NUMNP, NUMEL, NELBLK
      INTEGER NUMNPS, LNPSNL, LNPSDF
      INTEGER NUMESS, LESSEL, LESSNL, LESSDF
      LOGICAL ALL

      ALL = (OPTION .EQ. '*')
      

      IF ((ALL) .OR. (INDEX (OPTION, 'I') .GT. 0)) THEN
         WRITE(*, 10000)
      ENDIF

      IF ((ALL) .OR. (INDEX (OPTION, 'N') .GT. 0)) THEN
         INQUIRE (NDB, NAME=FILNAM)
         WRITE (*, 10010) FILNAM(:LENSTR(FILNAM))
      END IF

      IF ((ALL) .OR. (INDEX (OPTION, 'T') .GT. 0)) THEN
         WRITE (*, 10020) TITLE(1:LENSTR(TITLE))
      END IF

      IF ((ALL) .OR. (INDEX (OPTION, 'I') .GT. 0)) THEN
         WRITE (*, 10030, IOSTAT=IDUM)
     &      NDIM, NUMNP, NUMEL, NELBLK
      END IF

      IF ((ALL) .OR. (INDEX (OPTION, 'S') .GT. 0)) THEN
         If (NUMNPS .LE. 0) THEN
            WRITE (*, 10040, IOSTAT=IDUM) NUMNPS
         ELSE
            WRITE (*, 10040, IOSTAT=IDUM) NUMNPS, LNPSNL, LNPSDF
         END IF
         IF (NUMESS .LE. 0) THEN
            WRITE (*, 10050, IOSTAT=IDUM) NUMESS
         ELSE
            WRITE (*, 10050, IOSTAT=IDUM) NUMESS, LESSEL, LESSNL, LESSDF
         END IF
      END IF

      IF ((ALL) .OR. (INDEX (OPTION, 'V') .GT. 0)) THEN
         WRITE (*, 10060, IOSTAT=IDUM) NVARGL, NVARNP, NVAREL
      END IF

10000 FORMAT (/, 1X, 'DATABASE INITIAL VARIABLES')
10010 FORMAT (/, 1X, 'Database:  ', A)
10020 FORMAT (/, 1X, A)
10030 FORMAT (
     &      /, 1X, 'Number of coordinates per node       =', I10
     &      /, 1X, 'Number of nodes                      =', I10
     &      /, 1X, 'Number of elements                   =', I10
     &      /, 1X, 'Number of element blocks             =', I10
     &      )
10040 FORMAT (
     &      /, 1X, 'Number of node sets                  =', I10, :
     &      /, 1X, '   Length of node list               =', I10
     &      /, 1X, '   Length of distribution list       =', I10
     &      )
10050 FORMAT
     &      (  1X, 'Number of side sets                  =', I10, :
     &      /, 1X, '   Length of element list            =', I10
     &      /, 1X, '   Length of node list               =', I10
     &      /, 1X, '   Length of distribution list       =', I10
     &      )
10060 FORMAT (
     &      /, 1X, 'Number of global variables           =', I10
     &      /, 1X, 'Number of variables at each node     =', I10
     &      /, 1X, 'Number of variables at each element  =', I10
     &      )

      RETURN
      END
