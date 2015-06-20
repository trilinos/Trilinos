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
      SUBROUTINE DBIELB (NDB, OPTION, NELBS, NELBE,
     &   IDELB, NUMELB, NUMLNK, NUMATR,
     &   A, KLINK, KATRIB, *)
C=======================================================================
C$Id: dbielb.f,v 1.3 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: dbielb.f,v $
CRevision 1.3  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.2  1997/03/20 19:40:12  caforsy
CUpdated Imakefile for Imake 6.1.  Changed printing routines to handle
Clarger problems.
C
CRevision 1.1.1.1  1990/08/14 16:12:37  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:12:35  gdsjaar
c Initial revision
c 
c Revision 1.1  90/08/09  13:39:08  gdsjaar
c Initial revision
c 

C   --*** DBIELB *** (EXOLIB) Read database element blocks
C   --   Written by Amy Gilkey - revised 10/14/87
C   --
C   --DBIELB reads the element block information from the database.
C   --An error message is displayed if the end of file is read.
C   --
C   --The dynamic memory arrays LINK and ATRIB must be reserved
C   --if the connectivity and attributes are to be stored.  These arrays
C   --will be expanded by this routine to hold the new data.
C   --
C   --Parameters:
C   --   NDB - IN - the database file
C   --   OPTION - IN - ' ' to not store, '*' to store all, else store options:
C   --      'H' to store all header information
C   --      'I' to store block IDs
C   --      'C' to store connectivity (+ NUMLNK + NUMELB)
C   --      'A' to store attributes (+ NUMATR + NUMELB)
C   --   NELBS, NELBE - IN - the number of first and last element blocks
C   --      to read
C   --   IDELB - OUT - the element block IDs for each block (if OPTION)
C   --   NUMELB - OUT - the number of elements in each block (if OPTION)
C   --   NUMLNK - OUT - the number of nodes per element in each block
C   --      (if OPTION)
C   --   NUMATR - OUT - the number of attributes in each block (if OPTION)
C   --   A - IN/OUT - the dynamic memory base array
C   --   KLINK - OUT - pointer to the connectivity for each block (if OPTION)
C   --   KATRIB - OUT - pointer to the attributes for each block (if OPTION)
C   --   * - return statement if end of file or read error
C   --
C   --Database must be positioned at start of element block information
C   --upon entry; upon exit at end of element block information.

      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER NELBS, NELBE
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER NUMATR(*)
      DIMENSION A(*)
      INTEGER KLINK, KATRIB

      CHARACTER*80 ERRMSG

      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'C') .GT. 0)) THEN
         CALL MDFIND ('LINK', KLINK, IELNK)
      ELSE
         KLINK = 0
      END IF
      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'A') .GT. 0)) THEN
         CALL MDFIND ('ATRIB', KATRIB, IEATR)
      ELSE
         KATRIB = 0
      END IF
      CALL MDSTAT (NERR, MEM)
      IF (NERR .GT. 0) GOTO 110

      DO 100 NELB = NELBS, NELBE
         IELB = NELB-NELBS+1

         READ (NDB, END=120, ERR=120, IOSTAT=IERR)
     &      ID, NEL, NLNK, NATR
         IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'H') .GT. 0)) THEN
            IDELB(IELB) = ID
            NUMELB(IELB) = NEL
            NUMLNK(IELB) = NLNK
            NUMATR(IELB) = NATR
         END IF
         IF (INDEX (OPTION, 'I') .GT. 0) THEN
            IDELB(IELB) = ID
         END IF

         IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'C') .GT. 0)) THEN
            NUMELB(IELB) = NEL
            NUMLNK(IELB) = NLNK
            ISLNK = IELNK + 1
            IELNK = IELNK + NLNK * NEL
            CALL MDLONG ('LINK', KLINK, IELNK)
         END IF
         IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'A') .GT. 0)) THEN
            NUMELB(IELB) = NEL
            NUMATR(IELB) = NATR
            ISATR = IEATR + 1
            IEATR = IEATR + NATR * NEL
            CALL MDLONG ('ATRIB', KATRIB, IEATR)
         END IF
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) GOTO 110

         CALL DBIEB1 (NDB, OPTION, NELB,
     &      NEL, NLNK, NATR,
     &      A(KLINK+ISLNK-1), A(KATRIB+ISATR-1),
     $        MAX(NATR,1), MAX(NLNK,1), *130)
  100 CONTINUE

  110 CONTINUE
      RETURN

  120 CONTINUE
      WRITE (ERRMSG, '(A, I8)', IOSTAT=IDUM)
     &   'ELEMENT BLOCK SIZING PARAMETERS for BLOCK', NELB
      CALL DBERR (IERR, ERRMSG)
  130 CONTINUE
      RETURN 1
      END
