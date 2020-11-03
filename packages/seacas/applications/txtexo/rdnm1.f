C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE RDNM1 (NTXT, NDB, NELBLK, NVAREL, ISEVOK, LSEVOK, *)
C=======================================================================

C   --*** RDNM1 *** (TXTEXO) Internal to RDNAME
C   --   Written by Amy Gilkey - revised 02/22/88
C   --
C   --RDNM1 reads the element block variable truth table.
C   --
C   --Parameters:
C   --   NTXT - IN - the text file
C   --   NELBLK - IN - the number of element blocks
C   --   NVAREL - IN - the number of element variables
C   --   ISEVOK - OUT - the element block variable truth table;
C   --      variable i of block j exists iff ISEVOK(j,i)
C   --   * - return statement if error encountered, including end-of-file;
C   --      NO message is printed
C   --
C   --Database must be positioned in front of truth table upon entry;
C   --upon exit positioned after table.

      INTEGER ISEVOK(NVAREL,*)
      LOGICAL LSEVOK(NVAREL,*)

      CHARACTER*32 STRA

C ... Nothing to read if NVAREL == 0
      if (nvarel .eq. 0) return

      IELB = 0
      READ (NTXT, *, END=110, ERR=110)
      DO 100 IELB = 1, NELBLK
        READ (NTXT, *, END=110, ERR=110) (LSEVOK(I,IELB), I=1,NVAREL)
        DO 90 I = 1, NVAREL
          if (lsevok(i,ielb)) then
            isevok(i,ielb) = 1
          else
            isevok(i,ielb) = 0
          end if
 90     continue
 100  CONTINUE

      call expvtt(ndb, nelblk, nvarel, isevok, ierr)

      RETURN

 110  CONTINUE
      CALL INTSTR (1, 0, IELB, STRA, LSTRA)
      CALL PRTERR ('FATAL',
     &  'Reading ELEMENT BLOCK VARIABLE TRUTH TABLE for block '
     &  // STRA(:LSTRA))
      GOTO 120
 120  CONTINUE
      RETURN 1
      END
