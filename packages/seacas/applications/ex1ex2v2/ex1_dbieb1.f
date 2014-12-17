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
      SUBROUTINE DBIEB1 (NDB, OPTION, IELB, NUMELB, NUMLNK, NUMATR,
     &   LINK, ATRIB, NATRDM, NLNKDM, *)
C=======================================================================
C$Id: dbieb1.f,v 1.3 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: dbieb1.f,v $
CRevision 1.3  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.2  1997/03/20 19:40:11  caforsy
CUpdated Imakefile for Imake 6.1.  Changed printing routines to handle
Clarger problems.
C
CRevision 1.1.1.1  1990/08/14 16:12:34  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:12:33  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:08  gdsjaar
c Initial revision
c 

C   --*** DBIEB1 *** (EXOLIB) Read database element block misc.
C   --   Written by Amy Gilkey - revised 10/14/87
C   --   Modified by Greg Sjaardema - 8/8/90
C   --      ---Removed MAX from Dimension statements, Added NATRDM, NLNKDM
C   --
C   --DBIEB1 reads the element block connectivity and attribute information
C   --from the database.  An error message is displayed if the end of file
C   --is read.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   OPTION - IN - ' ' to not store, '*' to store all, else store options:
C   --      'C' to store connectivity
C   --      'A' to store attributes
C   --   IELB - IN - the element block number
C   --   NUMELB - IN - the number of elements in the block
C   --   NUMLNK - IN - the number of nodes per element;
C   --      negate to not store connectivity
C   --   NUMATR - IN - the number of attributes;
C   --      negate to not store attributes
C   --   LINK - OUT - the element connectivity for this block
C   --   ATRIB - OUT - the attributes for this block
C   --	 NATRDM - IN - dimension of atrib array
C   --   NLNKDM - IN - dimension of link array
C   --   * - return statement if end of file or read error
C   --
C   --Database must be positioned at start of element block misc. information
C   --upon entry; upon exit at end of element block misc. information.

      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER NUMELB, NUMLNK, NUMATR
      INTEGER LINK(NLNKDM, *)
      REAL ATRIB(NATRDM,*)

      CHARACTER*80 ERRMSG

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'C') .GT. 0)) THEN
         READ (NDB, END=100, ERR=100, IOSTAT=IERR)
     &      ((LINK(I,NE), I=1,NUMLNK), NE=1,NUMELB)
      ELSE
         READ (NDB, END=100, ERR=100, IOSTAT=IERR)
      END IF

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'A') .GT. 0)) THEN
         READ (NDB, END=110, ERR=110, IOSTAT=IERR)
     &      ((ATRIB(I,NE), I=1,NUMATR), NE=1,NUMELB)
      ELSE
         READ (NDB, END=110, ERR=110, IOSTAT=IERR)
      END IF

      RETURN

  100 CONTINUE
      WRITE (ERRMSG, '(A, I8)', IOSTAT=IDUM)
     &   'CONNECTIVITY for ELEMENT BLOCK', IELB
      GOTO 120
  110 CONTINUE
      WRITE (ERRMSG, '(A, I8)', IOSTAT=IDUM)
     &   'ATTRIBUTES for ELEMENT BLOCK', IELB
      GOTO 120
  120 CONTINUE
      CALL DBERR (IERR, ERRMSG)
      RETURN 1
      END
