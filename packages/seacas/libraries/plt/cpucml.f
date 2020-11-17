C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CPUCML(LINE,PROMPT,L)
      CHARACTER*(*) LINE,PROMPT
      LOGICAL FIRST
      DATA FIRST/.TRUE./

      IF (FIRST) THEN
         FIRST = .FALSE.
         LINE = ' '
         RETURN

      ELSE
         IF (L.LT.0) THEN
            CALL CHRTRM(PROMPT,LP)

         ELSE
            LP = L
         END IF

         IF (LP.NE.0) THEN
            WRITE (6,'(1x,a)') PROMPT(1:LP)
         END IF

         READ (5,'(a)',ERR=10,END=10) LINE
         CALL CHRTRM(LINE,L)
      END IF

      RETURN

   10 L = -1
      RETURN

      END
