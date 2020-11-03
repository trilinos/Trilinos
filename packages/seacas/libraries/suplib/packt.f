C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      SUBROUTINE PACKT (TITLE,LENGTH)

C ... REMOVE MULTIPLE BLANKS FROM A TITLE OR LABEL

      CHARACTER*(*) TITLE
      CHARACTER*1 BLANK
      DATA BLANK/' '/
      I=1
      L=1

C ... SKIP LEADING BLANKS

   10 CONTINUE
      IF (TITLE(I:I) .NE. BLANK) GO TO 20
      I=I+1
      GO TO 10
   20 CONTINUE
      TITLE(L:L)=TITLE(I:I)
   30 CONTINUE
      IF (I .GE. LENGTH) GO TO 60
      I=I+1
      L=L+1
      TITLE(L:L)=TITLE(I:I)
      IF (TITLE(I:I) .EQ. BLANK) THEN
   40    CONTINUE
         IF (TITLE(I:I) .NE. BLANK .OR. I .GE. LENGTH) GO TO 50
         I=I+1
         GO TO 40
   50    CONTINUE
         L=L+1
         TITLE(L:L)=TITLE(I:I)
       END IF
      GO TO 30
   60 CONTINUE
      DO 70 I=L+1,LENGTH
         TITLE(I:I)=BLANK
   70    CONTINUE
      END
