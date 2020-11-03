C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SELELM (MAT, SELECT, NUMEL, NELBLK, NUMSEL)
C=======================================================================
      DIMENSION MAT(6,*)
      LOGICAL SELECT(*)

      CALL INILOG (NUMEL, .FALSE., SELECT)

      DO 20 IBLK = 1, NELBLK
         IF (MAT(5,IBLK) .GT. 0) THEN
            IBEG = MAT(3,IBLK)
            IEND = MAT(4,IBLK)
            DO 10 IEL = IBEG, IEND
               SELECT(IEL) = .TRUE.
   10       CONTINUE
         END IF
   20 CONTINUE

      NUMSEL = NUMEQL (.TRUE., NUMEL, SELECT)

      RETURN
      END
