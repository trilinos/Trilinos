C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE LINEPR (ML, MP, LINKP, LCON, II, I1, I2, I3, J1, J2,
     &   J3)
C***********************************************************************

C  SUBROUTINE LINEPR = GETS THE LINE PARAMETERS
C***********************************************************************

      DIMENSION LCON(3, ML)
      DIMENSION LINKP(2, MP)

      LOGICAL ADDLNK

      ADDLNK = .FALSE.
      I1 = LCON (1, II)
      I2 = LCON (2, II)
      I3 = LCON (3, II)
      CALL LTSORT (MP, LINKP, I1, J1, ADDLNK)
      CALL LTSORT (MP, LINKP, I2, J2, ADDLNK)
      IF (I3 .NE. 0) THEN
         CALL LTSORT (MP, LINKP, IABS(I3), J3, ADDLNK)
      ELSE
         J3 = 0
      ENDIF

      RETURN

      END
