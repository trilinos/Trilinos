C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SELNOD (MAT, IX, SELECT, NUMNP, NNODES, NELBLK, NUMSEL)
C=======================================================================
      DIMENSION MAT(6,*), IX(NNODES,*)
      LOGICAL SELECT(*)

      CALL INILOG (NUMNP, .FALSE., SELECT)

      DO 30 IBLK = 1, NELBLK
         IF (MAT(5,IBLK) .GT. 0) THEN
            IBEG = MAT(3,IBLK)
            IEND = MAT(4,IBLK)
            DO 20 IEL = IBEG, IEND
               DO 10 INOD = 1, NNODES
                  SELECT(IX(INOD, IEL)) = .TRUE.
   10          CONTINUE
   20       CONTINUE
         END IF
   30 CONTINUE

      NUMSEL = NUMEQL (.TRUE., NUMNP, SELECT)

      RETURN
      END
