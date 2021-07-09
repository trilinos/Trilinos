C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE GRLWID
C=======================================================================

C   --*** GRLWID *** (GRPLIB) Set line widths (PLT)
C   --   Written by Amy Gilkey - revised 05/26/87
C   --
C   --GRLWID sets all controllable line widths.  The line width is a
C   --constant value (determined by the device) multiplied by the ratio
C   --of the existing value to the last constant value.  The first time
C   --this routine is called, all widths are set to the constant value.
C   --
C   --Common Variables:
C   --   Uses ICURDV, DEVOK, DEVCOD, NSNAP of /GRPCOM/

C   --Routines Called:
C   --   PLTGTG - (PLTLIB) Get graph parameters (see PLTSTG)
C   --   PLTGTT - (PLTLIB) Get text parameters (see PLTSTT)
C   --   PLTGTV - (PLTLIB) Get vector parameters (see PLTSTV)
C   --   PLTSTG - (PLTLIB) Set graph parameters:
C   --      28 = (KAXESZ) axis size
C   --      33 = (KTICSZ) tick-mark size
C   --      34 = (KGRISZ) grid-line size
C   --      29 = (KCRVSZ) line size
C   --   PLTSTT - (PLTLIB) Set text parameters:
C   --      11 = (KTXTSZ) text width
C   --   PLTSTV - (PLTLIB) Set vector parameters:
C   --      2 = (KLINSZ) line width

      PARAMETER (KAXESZ=28, KCRVSZ=29, KNUMSZ=30, KLABSZ=31,
     &   KSYMSZ=32, KTICSZ=33, KGRISZ=34, KMINSZ=35, KZLNSZ=36)
      PARAMETER (KTXTSZ=11)
      PARAMETER (KLINSZ=2)

      COMMON /GRPCOC/ DEVNAM(2), DEVCOD(2)
      CHARACTER*3 DEVNAM
      CHARACTER*8 DEVCOD
      COMMON /GRPCOM/ ICURDV, ISHARD, DEVOK(2), TALKOK(2),
     &   NSNAP(2), IFONT(2), SOFTCH(2), AUTOPL(2),
     &   MAXCOL(2), NUMCOL(0:1,2), MAPALT(2), MAPUSE(2)
      LOGICAL ISHARD, DEVOK, TALKOK, SOFTCH, AUTOPL

      LOGICAL LDUM, PLTGTT, PLTSTT, PLTGTG, PLTSTG, PLTGTV, PLTSTV

      LOGICAL FIRST
      SAVE FIRST

      REAL DEVSZ(2), FIXSIZ
      SAVE DEVSZ, FIXSIZ

      PARAMETER (NUMSTG = 9)
      INTEGER KSTG(NUMSTG)
      SAVE KSTG

      DATA FIRST / .TRUE. /
      DATA KSTG / KAXESZ, KCRVSZ, KNUMSZ, KLABSZ,
     &   KSYMSZ, KTICSZ, KGRISZ, KMINSZ, KZLNSZ /

      IF (FIRST) THEN
         DO 100 IDEV = 1, 2
            IF (DEVOK(IDEV)) THEN
               IF (DEVCOD(IDEV) .EQ. 'DICOMED') THEN
                  DEVSZ(IDEV) = 10.0
               ELSE
                  DEVSZ(IDEV) = 160.0
               END IF
            END IF
  100    CONTINUE

         OLDRAT = 1.0
         FIXSIZ = 0.0
      END IF

      OLDFIX = FIXSIZ
      FIXSIZ = DEVSZ(ICURDV)

      IF (.NOT. FIRST) THEN
         LDUM = PLTGTT (KTXTSZ, OLDSIZ)
         OLDRAT = OLDSIZ / OLDFIX
      END IF
      LDUM = PLTSTT (KTXTSZ, FIXSIZ*OLDRAT)
      DO 110 I = 1, NUMSTG
         IF (.NOT. FIRST) THEN
            LDUM = PLTGTG (KSTG(I), OLDSIZ)
            OLDRAT = OLDSIZ / OLDFIX
         END IF
         LDUM = PLTSTG (KSTG(I), FIXSIZ*OLDRAT)
  110 CONTINUE
      IF (.NOT. FIRST) THEN
         LDUM = PLTGTV (KLINSZ, OLDSIZ)
         OLDRAT = OLDSIZ / OLDFIX
      END IF
      LDUM = PLTSTV (KLINSZ, FIXSIZ*OLDRAT)

      FIRST = .FALSE.

      RETURN
      END
