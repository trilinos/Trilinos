C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE INTSTR (NNUM, LNGSTR, INUM, ISTR, LSTR)
C=======================================================================
C$Id: intstr.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: intstr.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:15:14  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:15:12  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:32  gdsjaar
c Initial revision
c

C   --*** INTSTR *** (STRLIB) Convert integer numbers to strings
C   --   Written by Amy Gilkey - revised 03/14/88
C   --
C   --INTSTR converts a set of integer numbers into a consistent set of
C   --strings of the same length (right justified).
C   --
C   --Parameters:
C   --   NNUM - IN - the number of integers numbers in the set
C   --   LNGSTR - IN - the number of digits in the string; <= 0 if length
C   --      of longest integer
C   --   INUM - IN - the array of integers to be converted
C   --   ISTR - OUT - the set of integer strings
C   --   LSTR - OUT - the maximum length of the number strings

      INTEGER NNUM
      INTEGER LNGSTR
      INTEGER INUM(*)
      CHARACTER*(*) ISTR(*)
      INTEGER LSTR

      CHARACTER*20 SCISTR
      CHARACTER*5 FFMT

      IF ((NNUM .EQ. 1) .AND. (LNGSTR .LE. 0)) THEN

C      --Handle special case of single number, smallest length

         WRITE (FFMT, 10000, IOSTAT=IDUM) LEN (ISTR(1))

         WRITE (ISTR(1), FFMT, IOSTAT=IDUM) INUM(1)
         CALL SQZSTR (ISTR(1), LSTR)

      ELSE

C      --Determine smallest length of largest number string

         IF (LNGSTR .LE. 0) THEN
            MINE = 0
            MAXE = 0
            DO 100 I = 1, NNUM
               MINE = MIN (MINE, INUM(I))
               MAXE = MAX (MAXE, INUM(I))
  100       CONTINUE

            MM = MAX (IABS (MAXE), IABS (MINE))
            IF (MM .EQ. IABS (MINE)) MM = -MM
            WRITE (SCISTR, '(I12)', IOSTAT=IDUM) MM
            CALL SQZSTR (SCISTR, LSTR)

         ELSE
            LSTR = LNGSTR
         END IF

         LSTR = MIN (LSTR, LEN (ISTR(1)))

C      --Convert all integers to string of given length

         WRITE (FFMT, 10000, IOSTAT=IDUM) LSTR
10000     FORMAT ('(I', I2.2, ')')

         DO 110 I = 1, NNUM
            WRITE (ISTR(I), FFMT, IOSTAT=IDUM) INUM(I)
  110    CONTINUE
      END IF

      RETURN
      END
