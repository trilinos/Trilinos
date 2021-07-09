C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GRSYMB (LINTYP, ISYTYP, INDX)
C=======================================================================

C   --*** GRSYMB *** (GRPLIB) Set line type and symbol (PLT)
C   --   Written by Amy Gilkey - revised 02/14/86
C   --
C   --GRSYMB must first be called with a negative index to initialize
C   --the number of colors to be used (including black and white).
C   --
C   --GRSYMB sets the line type and symbol for PLTGRH depending on the
C   --passed index.  The line type and symbol are chosen consecutively
C   --(see LINTYP, ISYTYP), wrapping around if necessary.
C   --
C   --Parameters:
C   --   LINTYP - IN - the line type:
C   --      0 = no line
C   --     -n = vary line type (on INDX)
C   --      n = line type n
C   --   ISYTYP - IN - the symbol type:
C   --      0 = no symbol
C   --     -n = vary symbol type (on INDX)
C   --      n = symbol n
C   --   INDX - IN - the line type and symbol index

C   --Routines Called:
C   --   PLTSTG - (PLTLIB) Set graph parameter
C   --      6 = (KCOLIN) set line color for PLTGRH lines
C   --      44 = (KCOSYM) set symbol color for PLTGRH lines
C   --      5 = (KONLIN) set line type for PLTGRH lines
C   --      7 = (KONSYM) set symbol for PLTGRH lines

      PARAMETER (KCOLOR=1)
      PARAMETER (KCOLIN=6, KCOSYM=44, KONLIN=5, KONSYM=7)

      PARAMETER (NUMLIN=6, NUMSYM=6)
C      --NUMLIN - the number of line types
C      --NUMSYM - the number of symbols

      INTEGER INDX

      LOGICAL LDUM, PLTSTG
      INTEGER LSTLIN, LSTSYM
      SAVE LSTLIN, LSTSYM
C      --LSTLIN - the last line type set
C      --LSTSYM - the last symbol set

      DATA LSTLIN, LSTSYM / -999, -999 /

      IF ((LINTYP .LE. -1) .OR. (LINTYP .NE. LSTLIN)) THEN
         IF (LINTYP .LE. -1) THEN
            ILIN = MOD (INDX-1, NUMLIN) + 1
         ELSE
            ILIN = MIN (LINTYP, NUMLIN)
         END IF
         LDUM = PLTSTG (KONLIN,1.0*ILIN)
         LSTLIN = ILIN
      END IF
      IF ((ISYTYP .LE. -1) .OR. (ISYTYP .NE. LSTSYM)) THEN
         IF (ISYTYP .LE. -1) THEN
            ISYM = MOD (INDX-1, NUMSYM) + 1
         ELSE
            ISYM = MIN (ISYTYP, NUMSYM)
         END IF
         LDUM = PLTSTG (KONSYM, 1.0*ISYM)
         LSTSYM = ISYM
      END IF

      RETURN
      END
