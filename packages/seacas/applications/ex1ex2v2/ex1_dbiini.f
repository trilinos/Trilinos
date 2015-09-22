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
      SUBROUTINE DBIINI (NDB, OPTION, NVERS, TITLE,
     &   NDIM, NUMNP, NUMEL, NELBLK,
     &   NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL, *)
C=======================================================================
C$Id: dbiini.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: dbiini.f,v $
CRevision 1.2  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:12:43  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:12:42  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:09  gdsjaar
c Initial revision
c 

C   --*** DBIINI *** (EXOLIB) Read database title and initial variables
C   --   Written by Amy Gilkey - revised 05/24/88
C   --
C   --DBIINI rewinds the database and reads the title and the initial
C   --variables from the database.  An error message is displayed if
C   --the end of file is read.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   OPTION - IN - ' ' to not store, '*' to store all, else store options:
C   --      'T' to store title
C   --      'I' to store initial variables
C   --   NVERS - OUT - the version number (not read, always 1)
C   --   TITLE - OUT - the database title (if OPTION)
C   --   NDIM - OUT - the number of coordinates per node
C   --   NUMNP - OUT - the number of nodes
C   --   NUMEL - OUT - the number of elements
C   --   NELBLK - OUT - the number of element blocks
C   --   NUMNPS - OUT - the number of node sets
C   --   LNPSNL - OUT - the length of the node sets node list
C   --   NUMESS - OUT - the number of side sets
C   --   LESSEL - OUT - the length of the side sets element list
C   --   LESSNL - OUT - the length of the side sets node list
C   --   * - return statement if end of file or read error
C   --
C   --Database is rewound upon entry; upon exit positioned at end of
C   --initial variables.

      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER NVERS
      CHARACTER*80 TITLE
      INTEGER NDIM, NUMNP, NUMEL, NELBLK,
     &   NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL

      CHARACTER*80 ERRMSG

      REWIND (NDB)

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'T') .GT. 0)) THEN
         READ (NDB, END=100, ERR=100, IOSTAT=IERR) TITLE
      ELSE
         READ (NDB, END=100, ERR=100, IOSTAT=IERR)
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'I') .GT. 0)) THEN
         READ (NDB, END=110, ERR=110, IOSTAT=IERR)
     &      NUMNP, NDIM, NUMEL, NELBLK,
     &      NUMNPS, LNPSNL, NUMESS, LESSEL, LESSNL
         NVERS = 1
      ELSE
         READ (NDB, END=110, ERR=110, IOSTAT=IERR)
      END IF

      RETURN

  100 CONTINUE
      ERRMSG = 'TITLE'
      GOTO 120
  110 CONTINUE
      ERRMSG = 'INITIAL VARIABLES'
      GOTO 120
  120 CONTINUE
      CALL DBERR (IERR, ERRMSG)
      RETURN 1
      END
