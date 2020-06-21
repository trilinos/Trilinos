C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    See packages/seacas/LICENSE for details

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
