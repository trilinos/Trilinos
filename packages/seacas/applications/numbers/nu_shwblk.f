C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE SHWBLK (NELBLK, MAT, NSELND, NSELEL)
      DIMENSION MAT(6, NELBLK)
      CHARACTER*16 TYPE
      CHARACTER*80 STRTMP
      include 'nu_io.blk'

      DO 10 IO=IOMIN, IOMAX
         WRITE (IO, 50)
   10 CONTINUE
      DO 30 ITMP=1, NELBLK
         I = MAT(6, ITMP)
         IF (MAT(5,I) .EQ. 1) THEN
            TYPE = 'Selected'
         ELSE
            TYPE = 'Not Selected'
         END IF
         DO 20 IO=IOMIN, IOMAX
            WRITE (IO, 60) MAT(1,I), MAT(3,I), MAT(4,I), MAT(2,I), I,
     *         TYPE
   20    CONTINUE
   30 CONTINUE
      WRITE (STRTMP, 70) NSELND, NSELEL
      CALL SQZSTR(STRTMP, LSTR)
      DO 40 IO=IOMIN, IOMAX
         WRITE (IO, 80) STRTMP(:LSTR)
   40 CONTINUE

      RETURN
   50 FORMAT (//
     *   '     Material  First     Last   Number of     Block'/,
     *   '      Number  Element   Element  Elements ')
   60 FORMAT (5I10,5X,A16)
   70 FORMAT (I8,' nodes and ',I8,' elements selected')
   80 FORMAT (/5X,A)
      END
