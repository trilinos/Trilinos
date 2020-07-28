C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE DBIELB (NDB, OPTION, NELBS, NELBE, IDELB, NUMELB,
     &           NUMLNK, NUMATR, A, IA, KLINK, KATRIB, NAMELB, *)
C=======================================================================
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
C   --   NAMELB - OUT - the type of element in each block
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
      CHARACTER*(MXSTLN) NAMELB(*)

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
         call exgelb(ndb, idelb(ielb), namelb(ielb), nel, nlnk,
     &               natr, ierr)
         IELNK = IELNK + NLNK * NEL
         IEATR = IEATR + NATR * NEL
 90    CONTINUE
       IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'C') .GT. 0)) THEN
         CALL MDLONG ('LINK', KLINK, IELNK)
       END IF
       IF ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'A') .GT. 0)) THEN
         CALL MDLONG ('ATRIB', KATRIB, IEATR)
       END IF
       CALL MDSTAT (NERR, MEM)
       IF (NERR .GT. 0) GOTO 110
      END IF

      IELNK = IESAV
      IEATR = IASAV
      DO 100 NELB = NELBS, NELBE
         IELB = NELB-NELBS+1
         call exgelb(ndb, idelb(ielb), namelb(ielb), nel, nlnk,
     &               natr, ierr)
         CALL EXUPCS (NAMELB(IELB))
         CALL PCKSTR (1, NAMELB(IELB))

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
      ISIZE = 0
      DO 105 I = NELBS, NELBE
         IELB = I-NELBS+1
         ISIZE = ISIZE + NUMELB(IELB)*NUMLNK(IELB)
  105 CONTINUE
  110 CONTINUE
      RETURN

  130 CONTINUE
      RETURN 1

      END
