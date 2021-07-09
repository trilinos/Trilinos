C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
C=======================================================================
      SUBROUTINE DBIELB (NDB, OPTION, NELBS, NELBE, IDELB, NUMELB,
     &           NUMLNK, NUMATR, BLKTYP, A, IELNK, IEATR, IOERR)
C=======================================================================

C   --*** DBIELB *** (EXOLIB) Read database element blocks
C   --   Written by Amy Gilkey - revised 10/14/87
C   --   Modified to Read EXODUSIIV2
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
C   --                  Element block ID's are always read in
C   --                  'H' to store all header information
C   --                  'C' to store connectivity (+ NUMLNK + NUMELB)
C   --                  'A' to store attributes (+ NUMATR + NUMELB)
C   --   NELBS  - IN  - the number of the first element block to read
C   --   NELBE  - IN  - the number of the last element block to read
C   --   IDELB  - OUT - the element block IDs for each block
C   --   NUMELB - OUT - the number of elements in each block (if OPTION)
C   --   NUMLNK - OUT - the number of nodes per element in each block
C   --                  (if OPTION)
C   --   NUMATR - OUT - the number of attributes in each block (if OPTION)
C   --   BLKTYP - OUT - the type of elements in the each block
C   --   A      - I/O - the dynamic memory base array
C   --   IELNK  - OUT - the size of the connectivity array
C   --   IEATR  - OUT - the size of the attribute array
C   --   IOERR  - OUT - error flag
C   --

      include 'exodusII.inc'

      INTEGER NDB
      CHARACTER*(*) OPTION
      INTEGER NELBS, NELBE
      INTEGER IDELB(*)
      INTEGER NUMELB(*)
      INTEGER NUMLNK(*)
      INTEGER NUMATR(*)
      DIMENSION A(1)
      INTEGER KLINK, KATRIB
      CHARACTER*(MXSTLN) BLKTYP(*)
      INTEGER IELNK, IEATR
      INTEGER IOERR
      LOGICAL ALL, HOPT, COPT, AOPT

      IOERR = 0
      IELNK = 0
      IEATR = 0

C     Set option flags
      ALL  = (OPTION .EQ. '*')
      HOPT = ALL .OR. (INDEX (OPTION, 'H') .GT. 0)
      COPT = ALL .OR. (INDEX (OPTION, 'C') .GT. 0)
      AOPT = ALL .OR. (INDEX (OPTION, 'A') .GT. 0)
C     No need to check if option = 'I' because the element block
C     IDs are always read in this subroutine.  The header, connectivity
C     and attribute arrays need the element block IDs.

C     Read the element block ID's
      CALL EXGEBI (NDB, IDELB, IERR)

      IF (ALL .OR. HOPT) THEN
         DO 100 IELB = NELBS, NELBE
C           Read element block parameters - returns:
C           1. element block ID's
C           2. element type in the element block
C           3. number of elements in the element block
C           4. number of nodes per element in the element block
C           5. number of attributes per element in the element block
C           6. error id
            call exgelb(ndb, idelb(ielb), blktyp(ielb), numelb(ielb),
     &                  numlnk(ielb), numatr(ielb), ierr)
            call exupcs(blktyp(ielb))
            call pckstr(1, blktyp(ielb))
  100    CONTINUE
      END IF

      IF (ALL .OR. COPT) THEN
C        Assumption: space has been reserved for 'LINK'
C        call MxFIND(array_name, ret_index, array_size)
         CALL MDFIND ('LINK', KLINK, IELNK)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) then
            CALL MEMERR
            IOERR = 1
            RETURN
         END IF
         DO 110 IELB = NELBS, NELBE
C           number of elements in a block
            NEL  = NUMELB(IELB)
C           number of nodes per element in the element block
            NLNK = NUMLNK(IELB)
            ISLNK = IELNK + 1
            IELNK = IELNK + NLNK * NEL
            CALL MDLONG ('LINK', KLINK, IELNK)
            CALL EXGELC(NDB, IDELB(IELB), A(KLINK+ISLNK-1), IERR)
  110    CONTINUE
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) then
            CALL MEMERR
            IOERR = 1
            RETURN
         END IF
      END IF

      IF (ALL .OR. AOPT) THEN
C        Assumption: space has been reserved for 'ATRIB'
C        call MxFIND(array_name, ret_index, array_size)
         CALL MDFIND ('ATRIB', KATRIB, IEATR)
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) then
            CALL MEMERR
            IOERR = 1
            RETURN
         END IF
         DO 120 IELB = NELBS, NELBE
C           number of elements in a block
            NEL  = NUMELB(IELB)
C           number of attributes in this block
            NATR = NUMATR(IELB)
            IF (NATR .GT. 0) THEN
               ISATR = IEATR + 1
               IEATR = IEATR + NATR * NEL
               CALL MDLONG ('ATRIB', KATRIB, IEATR)
               CALL EXGEAT(NDB, IDELB(IELB), A(KATRIB+ISATR-1), IERR)
            END IF
  120    CONTINUE
         CALL MDSTAT (NERR, MEM)
         IF (NERR .GT. 0) then
            CALL MEMERR
            IOERR = 1
         END IF
      END IF

      RETURN
      END
