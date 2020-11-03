C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE LISSEL (OPT, TYPE, IOMIN, IOMAX, LIST, SELECT, NUMLST)
C=======================================================================

C ... Output selected entities to Terminal and/or List file

C     OPT = IN = Option:
C                 'L' = Selected by logical list
C                 'A' = All selected are in list
C                 'R' = List in Range form

      LOGICAL SELECT(*)
      CHARACTER*(*) OPT, TYPE
      INTEGER LIST(*), ISCR(12)
      CHARACTER*80 STRING
      LOGICAL ISABRT

      WRITE (STRING, 10) TYPE(:LENSTR(TYPE))
   10 FORMAT ('List of Selected ',A)
      CALL SQZSTR(STRING,LSTR)
      DO 20 IO=IOMIN, IOMAX
         WRITE (IO, 30) STRING(:LSTR)
   20 CONTINUE
   30 FORMAT (/,1X,A/)

      II = 0
      IF (INDEX(OPT,'R') .GT. 0) THEN
         CALL RANGE(NUMLST, SELECT, IOMIN, IOMAX)
      ELSE
      DO 50 I = 1, NUMLST
         IF (SELECT(I)) THEN
            II = II + 1
            ISCR(II) = I
            IF (II .EQ. 12) THEN
               IF (ISABRT()) RETURN
               DO 40 IO = IOMIN, IOMAX
                  WRITE (IO, 70) (ISCR(ILST),ILST=1,II)
   40          CONTINUE
               II = 0
            END IF
         END IF
   50 CONTINUE
      IF (II .GT. 0) THEN
         DO 60 IO = IOMIN, IOMAX
            WRITE (IO, 70) (ISCR(ILST),ILST=1,II)
   60    CONTINUE
      END IF
   70 FORMAT ((1X,12I6))
      END IF
      RETURN
      END
