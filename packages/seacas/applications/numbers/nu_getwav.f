C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GETWAV (MAT, WAVE, NELBLK, LABEL)
      DIMENSION MAT(6,*), WAVE(*)
      DIMENSION IDUM(4), RV(4), KV(4)
      CHARACTER*16 LABEL(*), CV(4)
      CHARACTER*32 PRMPT

      WRITE (*, 20)
   20 FORMAT (/,'   Input WAVESPEED, and NAME: ')
      DO 40 IBLK=1,NELBLK
         I = MAT(6,IBLK)
         WRITE (PRMPT, 30) MAT(1,I)
   30    FORMAT ('    Material ',I5,' > ')
         CALL FREFLD (0, 0, PRMPT(:LENSTR(PRMPT)+1), 2, IOS,
     *      NF, KV, CV, IDUM, RV)
         WAVE(I)    = RV(1)
         IF (KV(2) .NE. -1) LABEL(I) = CV(2)
   40 CONTINUE
      RETURN
      END
