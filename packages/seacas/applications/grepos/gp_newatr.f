C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE NEWATR (NELBLK, NUMATR, ATRSC, NUMELB, ATRIB)
C=======================================================================

      INTEGER NUMATR(*)
      INTEGER NUMELB(*)
      REAL ATRSC(2,*)
      REAL ATRIB(*)

      IEATR = 0
      IAT   = 1

      DO 100 IELB = 1, NELBLK
         ISATR = IEATR + 1
         IEATR = IEATR + NUMATR(IELB) * NUMELB(IELB)
         if (numatr(ielb) .gt. 0) then
           CALL NEWAT1 (NUMELB(IELB), NUMATR(IELB),
     *       ATRSC(1,IAT), ATRIB(ISATR))
           IAT = IAT + NUMATR(IELB)
         end if
  100 CONTINUE

      RETURN
      END

      SUBROUTINE NEWAT1(NUMEL, NUMATR, ATRSC, ATRIB)
      REAL ATRSC(2,*)
      REAL ATRIB(*)

      IBEG = 1
      DO 110 IATR = 1, NUMATR
        if (ATRSC(1,IATR) .NE. 0.0) then
          DO 100 IEL = 1, NUMEL
            ATRIB(IBEG+NUMATR*(IEL-1)) = ATRSC(1,IATR)
 100      CONTINUE
        else if (ATRSC(2,IATR) .NE. 1.0) then
          DO 105 IEL = 1, NUMEL
            ATRIB(IBEG+NUMATR*(IEL-1)) = ATRSC(2,IATR) *
     *        ATRIB(IBEG+NUMATR*(IEL-1))
 105      CONTINUE
        end if
        IBEG = IBEG + 1
 110  CONTINUE

      RETURN
      END

