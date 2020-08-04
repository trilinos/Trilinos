C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION CPUNAL(IU)
      CHARACTER*40 CPUHLB
      COMMON /CPHLBN/CPUHLB
      INTEGER CPUNIT(20)
      INTEGER CPUNIF
      COMMON /CPUN/CPUNIT,CPUNIF
      CHARACTER*20 TERM
      CHARACTER*80 LINBUF(24)
      COMMON /REBUF/LINBUF,TERM
      INTEGER LINLEN(24)
      COMMON /REBUF2/NL,LINLEN

      CPUNAL = .TRUE.
      IF (CPUNIF.NE.12345) THEN
         I = 1
 2260    IF (.NOT. (I.LE.20)) GO TO 2280
         CPUNIT(I) = 80 + I - 1
         I = I + 1
         GO TO 2260

 2280    CONTINUE
         CPUNIF = 12345
      END IF

      I = 1
 2290 IF (.NOT. (I.LE.20)) GO TO 2310
      IF (CPUNIT(I).EQ.0) THEN
         GO TO 2300

      END IF

      IU = CPUNIT(I)
      CPUNIT(I) = 0
      RETURN

 2300 I = I + 1
      GO TO 2290

 2310 CONTINUE
      CPUNAL = .FALSE.
      CALL CPUERR('Cannot allocate logical unit.',2)
      RETURN

      END
