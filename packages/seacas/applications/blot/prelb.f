C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE PRELB (OPTION, NOUT, NELBLK, NLISEL, LISEL,
     &   IDELB, LENE, NUMLNK, NUMATR, LINK, ATRIB,
     &   NAMELB, EBNAME, NVAREL, NAMEEV, ISEVOK, LISEV,
     *   MAPEL, MAPND)
C=======================================================================

C   --*** PRELB *** (BLOT) Display database element blocks
C   --   Written by Amy Gilkey - revised 01/18/88
C   --
C   --PRELB displays the element blocks.
C   --
C   --Parameters:
C   --   OPTION - IN - '*' to print all, else print options:
C   --      ' ' to print block summary
C   --      'N' to print element block name
C   --      'V' to print element variable truth table
C   --      'C' to print connectivity
C   --      'A' to print attributes
C   --   NOUT - IN - the output file, <=0 for standard
C   --   NELBLK - IN - the number of element blocks
C   --   NLISEL - IN - the number of selected elements by block
C   --   LISEL - IN - the indices of the selected elements by block
C   --   IDELB - IN - the element block ID for each block
C   --   LENE - IN - the cumulative element counts by element block
C   --   NUMLNK - IN - the number of nodes per element for each block
C   --   NUMATR - IN - the number of attributes for each block
C   --   LINK - IN - the connectivity array for all blocks
C   --   ATRIB - IN - the attribute array for all blocks
C   --   NAMELB - IN - the names of the element block types
C   --   NVAREL - IN - the number of element variables
C   --   NAMEEV - IN - the names of the element variables
C   --   ISEVOK - IN - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   LISEV - SCRATCH - size = NVAREL (if 'V' in OPTION)

      CHARACTER*(*) OPTION
      INTEGER NLISEL(0:*)
      INTEGER LISEL(0:*)
      INTEGER IDELB(*)
      INTEGER LENE(0:*)
      INTEGER NUMLNK(*)
      INTEGER NUMATR(*)
      INTEGER LINK(*)
      REAL ATRIB(*)
      CHARACTER*(*) NAMELB(*)
      CHARACTER*(*) EBNAME(*)
      CHARACTER*(*) NAMEEV(*)
      LOGICAL ISEVOK(NELBLK,NVAREL)
      INTEGER LISEV(*)
      INTEGER MAPEL(*), MAPND(*)

      LOGICAL ISABRT
      LOGICAL DONAM, DOVTBL, DOCONN, DOATR
      LOGICAL BLK1
      CHARACTER*32 STRA, STRB

      DONAM  = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'N') .GT. 0))
      DOVTBL = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'V') .GT. 0))
      DOCONN = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'C') .GT. 0))
      DOATR  = ((OPTION .EQ. '*') .OR. (INDEX (OPTION, 'A') .GT. 0))

      BLK1 = .TRUE.
      IF (NOUT .GT. 0) THEN
         IF (DOCONN .AND. DOATR) THEN
            WRITE (NOUT, 10020) 'CONNECTIVITY and ATTRIBUTES'
         ELSE IF (DOCONN) THEN
            WRITE (NOUT, 10020) 'CONNECTIVITY'
         ELSE IF (DOATR) THEN
            WRITE (NOUT, 10020) 'ATTRIBUTES'
         ELSE
            WRITE (NOUT, 10020)
         END IF
      END IF

      WRITE (STRA, 10000, IOSTAT=IDUM) NELBLK
10000  FORMAT ('(#', I5, ')')
      CALL PCKSTR (1, STRA)
      LSTRA = LENSTR (STRA)
      WRITE (STRB, 10010, IOSTAT=IDUM) LENE(NELBLK), LENE(NELBLK)
10010  FORMAT ('(', I10, '..', I10, ')')
      CALL PCKSTR (1, STRB)
      LSTRB = LENSTR (STRB)

      DO 110 IELB = 1, NELBLK
         IF (ISABRT ()) RETURN
         IF (NLISEL(IELB) .GT. 0) THEN
            IEL = LENE(IELB-1)+1
            LEL = LENE(IELB)

            IF (BLK1 .OR. DOCONN .OR. DOATR) THEN
               BLK1 = .FALSE.
               IF (NOUT .GT. 0) THEN
                  WRITE (NOUT, *)
               ELSE
                  WRITE (*, *)
               END IF
            END IF

            NUME = LENE(IELB) - LENE(IELB-1)
            WRITE (STRA, 10000, IOSTAT=IDUM) IELB
            CALL PCKSTR (1, STRA)
            WRITE (STRB, 10010, IOSTAT=IDUM) IEL, LEL
            CALL PCKSTR (1, STRB)
            IF (NOUT .GT. 0) THEN
               WRITE (NOUT, 10030, IOSTAT=IDUM)
     &            IDELB(IELB), STRA(:LSTRA),
     &            NUME, STRB(:LSTRB), NUMLNK(IELB), NUMATR(IELB)
            ELSE
               WRITE (*, 10030, IOSTAT=IDUM)
     &            IDELB(IELB), STRA(:LSTRA),
     &            NUME, STRB(:LSTRB), NUMLNK(IELB), NUMATR(IELB)
            END IF

            IF (DONAM) THEN
               LNAM = LENSTR(NAMELB(IELB))
               LNM =  LENSTR(EBNAME(IELB))
               IF (NOUT .GT. 0) THEN
                  WRITE (NOUT, 10040) EBNAME(IELB)(:LNM),
     $                 NAMELB(IELB)(:LNAM)
               ELSE
                  WRITE (*, 10040) EBNAME(IELB)(:LNM),
     $                 NAMELB(IELB)(:LNAM)
               END IF
            END IF

            IF (DOVTBL) THEN
               lname = 0
               NSEL = 0
               DO 100 I = 1, NVAREL
                  IF (ISEVOK(IELB,I)) THEN
                     l = lenstr(nameev(i))
                     if (l .gt. lname) lname = l
                     NSEL = NSEL + 1
                     LISEV(NSEL) = I
                  END IF
  100          CONTINUE

               IF (NOUT .GT. 0) THEN
                  WRITE (NOUT, 10050, IOSTAT=IDUM)
     &               (NAMEEV(LISEV(I))(:lname), I=1,NSEL)
               ELSE
                  WRITE (*, 10050, IOSTAT=IDUM)
     &               (NAMEEV(LISEV(I))(:lname), I=1,NSEL)
               END IF
            END IF

            IF (DOCONN .OR. DOATR) THEN
               IF (DOCONN) ISLNK = IDBLNK (IELB, 0, LENE, NUMLNK)
               IF (DOATR)  ISATR = IDBLNK (IELB, 0, LENE, NUMATR)
               CALL PREB1 (OPTION, NOUT, IEL-1,
     &            NLISEL(IELB), LISEL(IEL), NUMLNK(IELB), NUMATR(IELB),
     &            LINK(ISLNK), ATRIB(ISATR), MAPEL, MAPND)
               IF (ISABRT ()) RETURN
            END IF
         END IF
  110 CONTINUE

      RETURN

10020  FORMAT (/, 1X, 'ELEMENT BLOCKS', :, ' - ', A)
10030  FORMAT (1X, 'Block', I9, 1X, A, ':',
     &   I9, ' elements', 1X, A,
     &   I4, '-node', I4, ' attributes')
10040  FORMAT (4X, 'Element block name = "',A
     $      ,'", type = "', A, '"')
10050  FORMAT (4X, 'Defined variables:', :, 3X, 4 (2X, A), :, /
     &   (4X, 1X, 6 (2X, A)))
      END
