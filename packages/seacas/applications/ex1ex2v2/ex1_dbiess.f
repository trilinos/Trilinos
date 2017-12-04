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

C=======================================================================
      SUBROUTINE DBIESS (NDB, OPTION, NUMESS, LESSEL, LESSNL,
     &   IDESS, NEESS, NNESS, IXEESS, IXNESS, LTEESS, LTNESS, FACESS, *)
C=======================================================================
C$Id: dbiess.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: dbiess.f,v $
CRevision 1.2  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:12:40  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:12:38  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:09  gdsjaar
c Initial revision
c 

C   --*** DBIESS *** (EXOLIB) Read database node sets
C   --   Written by Amy Gilkey - revised 10/14/87
C   --
C   --DBIESS reads the side set information from the database.
C   --An error message is displayed if the end of file is read.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   OPTION - IN - ' ' to not store, '*' to store all, else store options:
C   --      'H' to store information about side sets
C   --      'E' to store side set elements
C   --      'N' to store side set nodes
C   --      'F' to store side set factors
C   --   NUMESS - IN - the number of side sets
C   --   LESSEL - IN - the length of the side sets element list
C   --   LESSNL - IN - the length of the side sets node list
C   --   IDESS - OUT - the side set ID for each set (if OPTION)
C   --   NEESS - OUT - the number of elements for each set (if OPTION)
C   --   NNESS - OUT - the number of nodes for each set (if OPTION)
C   --   IXEESS - OUT - the index of the first element for each set (if OPTION)
C   --   IXNESS - OUT - the index of the first node for each set (if OPTION)
C   --   LTEESS - OUT - the elements for all sets (if OPTION)
C   --   LTNESS - OUT - the nodes for all sets (if OPTION)
C   --   FACESS - OUT - the distribution factors for all sets (if OPTION)
C   --   * - return statement if end of file or read error
C   --
C   --Database must be positioned at start of side set information
C   --upon entry; upon exit at end of side set information.

      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER NUMESS, LESSEL, LESSNL
      INTEGER IDESS(*)
      INTEGER NEESS(*)
      INTEGER NNESS(*)
      INTEGER IXEESS(*)
      INTEGER IXNESS(*)
      INTEGER LTEESS(*)
      INTEGER LTNESS(*)
      REAL FACESS(*)

      CHARACTER*80 ERRMSG

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'H') .GT. 0)) THEN
         READ (NDB, END=100, ERR=100, IOSTAT=IERR)
     &      (IDESS(IESS), IESS=1,NUMESS)
         READ (NDB, END=110, ERR=110, IOSTAT=IERR)
     &      (NEESS(IESS), IESS=1,NUMESS)
         READ (NDB, END=120, ERR=120, IOSTAT=IERR)
     &      (NNESS(IESS), IESS=1,NUMESS)
         READ (NDB, END=130, ERR=130, IOSTAT=IERR)
     &      (IXEESS(IESS), IESS=1,NUMESS)
         READ (NDB, END=140, ERR=140, IOSTAT=IERR)
     &      (IXNESS(IESS), IESS=1,NUMESS)
      ELSE
         READ (NDB, END=100, ERR=100, IOSTAT=IERR)
         READ (NDB, END=110, ERR=110, IOSTAT=IERR)
         READ (NDB, END=120, ERR=120, IOSTAT=IERR)
         READ (NDB, END=130, ERR=130, IOSTAT=IERR)
         READ (NDB, END=140, ERR=140, IOSTAT=IERR)
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'E') .GT. 0)) THEN
         READ (NDB, END=150, ERR=150, IOSTAT=IERR)
     &      (LTEESS(NL), NL=1,LESSEL)
      ELSE
         READ (NDB, END=150, ERR=150, IOSTAT=IERR)
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'N') .GT. 0)) THEN
         READ (NDB, END=160, ERR=160, IOSTAT=IERR)
     &      (LTNESS(NL), NL=1,LESSNL)
      ELSE
         READ (NDB, END=160, ERR=160, IOSTAT=IERR)
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'F') .GT. 0)) THEN
         READ (NDB, END=170, ERR=170, IOSTAT=IERR)
     &      (FACESS(NL), NL=1,LESSNL)
      ELSE
         READ (NDB, END=170, ERR=170, IOSTAT=IERR)
      END IF

      RETURN

  100 CONTINUE
      ERRMSG = 'SIDE SET IDS'
      GOTO 180
  110 CONTINUE
      ERRMSG = 'SIDE SET NUMBER OF ELEMENTS'
      GOTO 180
  120 CONTINUE
      ERRMSG = 'SIDE SET NUMBER OF NODES'
      GOTO 180
  130 CONTINUE
      ERRMSG = 'SIDE SET ELEMENT INDICES'
      GOTO 180
  140 CONTINUE
      ERRMSG = 'SIDE SET NODE INDICES'
      GOTO 180
  150 CONTINUE
      ERRMSG = 'SIDE SET ELEMENTS'
      GOTO 180
  160 CONTINUE
      ERRMSG = 'SIDE SET NODES'
      GOTO 180
  170 CONTINUE
      ERRMSG = 'SIDE SET DISTRIBUTION FACTORS'
      GOTO 180
  180 CONTINUE
      CALL DBERR (IERR, ERRMSG)
      RETURN 1
      END
