C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CKCNTR (OK)
C=======================================================================

C   --*** CKCNTR *** (DETOUR) Check the contour values
C   --   Written by Amy Gilkey - revised 07/10/87
C   --
C   --CKCNTR checks that all specified contour values increase or decrease.
C   --
C   --Parameters:
C   --   OK - OUT - true iff the contour values are consistent
C   --
C   --Common Variables:
C   --   Uses CINTOK, LINCON, NCNTR, CINTV of /CNTR/

      include 'cntr.blk'
C     FLAG FOR EXACT CONTOUR VALUES FOR EACH PLOT
C     COMMON /CNTR/   CINTOK, LINCON, NCNTR, CMIN, CMAX, DELC,
C    &   CINTV(256), NOCMIN, NOCMAX, LABINC, MAXMIN, MAXMAX
C     LOGICAL CINTOK, LINCON, NOCMIN, NOCMAX

      LOGICAL OK

      CHARACTER*80 ERRSTR

      OK = .TRUE.

      IF (CINTOK) THEN
         IF (LINCON) THEN
            NC = NCNTR
         ELSE
            NC = NCNTR+1
         END IF

         IF (CINTV(1) .LE. CINTV(2)) THEN
            DO 100 I = 2, NC
               IF (CINTV(I-1) .GE. CINTV(I)) THEN
                  WRITE (ERRSTR, 10000)
     &               'Contour interval ', I-1, ' >= interval ', I
10000             FORMAT (A, I5, A, I5)
                  CALL SQZSTR (ERRSTR, LSTR)
                  CALL PRTERR ('CMDWARN', ERRSTR(:LSTR))
                  OK = .FALSE.
               END IF
  100       CONTINUE
         ELSE
            DO 110 I = 2, NC
               IF (CINTV(I-1) .LE. CINTV(I)) THEN
                  WRITE (ERRSTR, 10000)
     &               'Contour interval ', I-1, ' <= interval ', I
                  CALL SQZSTR (ERRSTR, LSTR)
                  CALL PRTERR ('CMDWARN', ERRSTR(:LSTR))
                  OK = .FALSE.
               END IF
  110       CONTINUE
         END IF
      END IF

      RETURN
      END
