C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SCLATR (NELBLK, IDELB, NUMATR, IDATR, DOALLA,
     *  IDBLK, DOALLB, ASCALE, ATRSCL, VERBOS, SCLSET)
C=======================================================================

      INTEGER IDELB(*)
      INTEGER NUMATR(*)
      LOGICAL DOALLA
      LOGICAL DOALLB
      REAL    ATRSCL(2,*)
C ... Row 1 == value to set attribute to, Row 2 == value to scale by.
      LOGICAL VERBOS
      LOGICAL SCLSET
C ... SCLSET == .TRUE. then do scale, .FALSE. then set to value

      IAT = 0
      DO 110 IELB = 1, NELBLK
        IF (IDELB(IELB) .EQ. IDBLK .OR. DOALLB) THEN
          DO 100 IATR = 1, NUMATR(IELB)
            IF (IDATR .EQ. IATR .OR. DOALLA) THEN
              IF (SCLSET) THEN
                ATRSCL(1,IAT+IATR) = 0.0
                ATRSCL(2,IAT+IATR) = ASCALE
                if (verbos) then
                  WRITE (*, 10000) IATR, IDELB(IELB), ' Scaled by ',
     *              ASCALE
                end if
              ELSE
                ATRSCL(1,IAT+IATR) = ASCALE
                ATRSCL(2,IAT+IATR) = 1.0
                if (verbos) then
                  WRITE (*, 10000) IATR, IDELB(IELB), ' Set to ', ASCALE
                end if
              END IF

            END IF
 100      CONTINUE
        END IF
        IAT = IAT + NUMATR(IELB)
 110  CONTINUE
      RETURN

10000 FORMAT (1X, 'Attribute', I3, ' of Block', I5,A, 1PE10.3)
      END

      SUBROUTINE CNTATR (NELBLK, NUMATR, IAT)
      INTEGER NUMATR(*)
      IAT = 0
      DO 110 IELB = 1, NELBLK
        IAT = IAT + NUMATR(IELB)
 110  CONTINUE
      RETURN
      END
