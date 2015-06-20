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
      SUBROUTINE DBINPS (NDB, OPTION, NUMNPS, LNPSNL,
     &   IDNPS, NNNPS, IXNNPS, LTNNPS, FACNPS, *)
C=======================================================================
C$Id: dbinps.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: dbinps.f,v $
CRevision 1.2  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:12:54  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:12:52  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:11  gdsjaar
c Initial revision
c 

C   --*** DBINPS *** (EXOLIB) Read database node sets
C   --   Written by Amy Gilkey - revised 10/14/87
C   --
C   --DBINPS reads the node set information from the database.
C   --An error message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   OPTION - IN - ' ' to not store, '*' to store all, else store options:
C   --      'H' to store information about node sets
C   --      'N' to store node set nodes
C   --      'F' to store node set factors
C   --   NUMNPS - IN - the number of node sets
C   --   LNPSNL - IN - the length of the node sets node list
C   --   IDNPS - OUT - the node set ID for each set (if OPTION)
C   --   NNNPS - OUT - the number of nodes for each set (if OPTION)
C   --   IXNNPS - OUT - the index of the first node for each set (if OPTION)
C   --   LTNNPS - OUT - the nodes for all sets (if OPTION)
C   --   FACNPS - OUT - the distribution factors for all sets (if OPTION)
C   --   * - return statement if end of file or read error
C   --
C   --Database must be positioned at start of node set information
C   --upon entry; upon exit at end of node set information.

      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER NUMNPS, LNPSNL
      INTEGER IDNPS(*)
      INTEGER NNNPS(*)
      INTEGER IXNNPS(*)
      INTEGER LTNNPS(*)
      REAL FACNPS(*)

      CHARACTER*80 ERRMSG

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'H') .GT. 0)) THEN
         READ (NDB, END=100, ERR=100, IOSTAT=IERR)
     &      (IDNPS(INPS), INPS=1,NUMNPS)
         READ (NDB, END=110, ERR=110, IOSTAT=IERR)
     &      (NNNPS(INPS), INPS=1,NUMNPS)
         READ (NDB, END=120, ERR=120, IOSTAT=IERR)
     &      (IXNNPS(INPS), INPS=1,NUMNPS)
      ELSE
         READ (NDB, END=100, ERR=100, IOSTAT=IERR)
         READ (NDB, END=110, ERR=110, IOSTAT=IERR)
         READ (NDB, END=120, ERR=120, IOSTAT=IERR)
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'N') .GT. 0)) THEN
         READ (NDB, END=130, ERR=130, IOSTAT=IERR)
     &      (LTNNPS(NL), NL=1,LNPSNL)
      ELSE
         READ (NDB, END=130, ERR=130, IOSTAT=IERR)
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'F') .GT. 0)) THEN
         READ (NDB, END=140, ERR=140, IOSTAT=IERR)
     &      (FACNPS(NL), NL=1,LNPSNL)
      ELSE
         READ (NDB, END=140, ERR=140, IOSTAT=IERR)
      END IF

      RETURN

  100 CONTINUE
      ERRMSG = 'NODE SET IDS'
      GOTO 150
  110 CONTINUE
      ERRMSG = 'NODE SET NUMBER OF NODES'
      GOTO 150
  120 CONTINUE
      ERRMSG = 'NODE SET INDICES'
      GOTO 150
  130 CONTINUE
      ERRMSG = 'NODE SET NODES'
      GOTO 150
  140 CONTINUE
      ERRMSG = 'NODE SET DISTRIBUTION FACTORS'
      GOTO 150
  150 CONTINUE
      CALL DBERR (IERR, ERRMSG)
      RETURN 1
      END
