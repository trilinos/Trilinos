C Copyright(C) 2011 Sandia Corporation.  Under the terms of Contract
C DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
C certain rights in this software
C 
C Redistribution and use in source and binary forms, with or without
C modification, are permitted provided that the following conditions are
C met:
C 
C * Redistributions of source code must retain the above copyright
C    notice, this list of conditions and the following disclaimer.
C           
C * Redistributions in binary form must reproduce the above
C   copyright notice, this list of conditions and the following
C   disclaimer in the documentation and/or other materials provided
C   with the distribution.
C                         
C * Neither the name of Sandia Corporation nor the names of its
C   contributors may be used to endorse or promote products derived
C   from this software without specific prior written permission.
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
      SUBROUTINE DBIELB (NDB, OPTION, NELBS, NELBE, IDELB, NUMELB,
     &           NUMLNK, NUMATR, A, IA, KLINK, KATRIB, BLKTYP, *)
C=======================================================================
C   --*** DBIELB *** (EXOLIB) Read database element blocks
C   --
C   --DBIELB reads the element block information from the database.
C   --An error message is displayed if the end of file is read.
C   --
C   --The dynamic memory arrays LINK and ATRIB must be reserved
C   --if the connectivity and attributes are to be stored.  These arrays
C   --will be expanded by this routine to hold the new data.
C   --
C   --Parameters:
C   --   NDB    - IN  - the database file
C   --   OPTION - IN  - ' ' to not store, '*' to store all, else store options:
C   --                  'H' to store all header information
C   --                  'I' to store block IDs
C   --                  'C' to store connectivity (+ NUMLNK + NUMELB)
C   --                  'A' to store attributes (+ NUMATR + NUMELB)
C   --   NELBS  - IN  - the number of the first element block to read
C   --   NELBE  - IN  - the number of the last element block to read
C   --   IDELB  - OUT - the element block IDs for each block 
C   --   NUMELB - OUT - the number of elements in each block 
C   --   NUMLNK - OUT - the number of nodes per element in each block
C   --   NUMATR - OUT - the number of attributes in each block 
C   --   A      - I/O - the dynamic memory base array for REALS
C   --   IA     - I/O - the dynamic memory base array for INTEGERS
C   --   KLINK  - OUT - pointer to the connectivity for each block 
C   --   KATRIB - OUT - pointer to the attributes for each block 
C   --   BLKTYP - OUT - the type of element in each block
C   --   *      - OUT - return statement if end of file or read error

      include 'exodusII.inc'
      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER NELBS, NELBE
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER NUMATR(*)
      DIMENSION A(*), IA(*)
      INTEGER KLINK, KATRIB
      CHARACTER*(MXSTLN) BLKTYP(*)

C ... Get element block ids
      call exgebi (ndb, idelb, ierr)

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

      IESAV = IELNK
      IASAV = IEATR
      IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'C') .GT. 0)
     &     .OR. (INDEX (OPTION, 'A') .GT. 0)) THEN
      DO 90 NELB = NELBS, NELBE
         IELB = NELB-NELBS+1
         call exgelb(ndb, idelb(ielb), blktyp(ielb), nel, nlnk,
     &               natr, ierr)
         IELNK = IELNK + NLNK * NEL
         IEATR = IEATR + NATR * NEL
 90    CONTINUE
       CALL MDLONG ('LINK', KLINK, IELNK)
       CALL MDLONG ('ATRIB', KATRIB, IEATR)
       CALL MDSTAT (NERR, MEM)
       IF (NERR .GT. 0) GOTO 110
      END IF

        
      IELNK = IESAV
      IEATR = IASAV 
      ISATR = 0
      ISLNK = 0
      DO 100 NELB = NELBS, NELBE
         IELB = NELB-NELBS+1
         call exgelb(ndb, idelb(ielb), blktyp(ielb), nel, nlnk,
     &               natr, ierr)
         CALL EXUPCS (BLKTYP(IELB))
         CALL PCKSTR (1, BLKTYP(IELB))

         IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'H') .GT. 0)) THEN
            NUMELB(IELB) = NEL
            NUMLNK(IELB) = NLNK
            NUMATR(IELB) = NATR
         END IF

         IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'C') .GT. 0)) THEN
            NUMELB(IELB) = NEL
            NUMLNK(IELB) = NLNK
            ISLNK = IELNK + 1
            IELNK = IELNK + NLNK * NEL
         END IF
         IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'A') .GT. 0)) THEN
            NUMELB(IELB) = NEL
            NUMATR(IELB) = NATR
            ISATR = IEATR + 1
            IEATR = IEATR + NATR * NEL
         END IF

         CALL DBIEBI (NDB, OPTION, IDELB(IELB), NEL, NLNK, NATR,
     &                IA(KLINK+ISLNK-1), A(KATRIB+ISATR-1),
     &                MAX(NATR,1), MAX(NLNK,1), *130)
  100 CONTINUE

C     Store the first pointers for each element block link array
      INC = 0
      ISIZE = 0
      DO 105 I = NELBS, NELBE
         IELB = I-NELBS+1
         INC = ISIZE + 1
         ISIZE = ISIZE + NUMELB(IELB)*NUMLNK(IELB)
  105 CONTINUE   
  110 CONTINUE
      RETURN

  130 CONTINUE
      RETURN 1

      END

