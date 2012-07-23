C
      LOGICAL FUNCTION CPUBRK(INPUT)
C***********************************************************************
C
C     FUNCTION CPUBRK = .TRUE. IF A CONTROL C HAS BEEN ENTERED AT TERMINAL
C
C***********************************************************************
      LOGICAL CPUIFC, INPUT
C
      IF (CPUIFC(INPUT)) THEN
         CPUBRK = .TRUE.
      ELSE
         CPUBRK = .FALSE.
      ENDIF
C
      RETURN
C
      END
