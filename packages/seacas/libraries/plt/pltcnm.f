C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      LOGICAL FUNCTION PLTCNM(VALUE,CLR)
      CHARACTER*(*) CLR

      PLTCNM = .TRUE.
      IF (VALUE.EQ.0.) THEN
         CLR = 'BLACK'

      ELSE IF (VALUE.EQ.1.) THEN
         CLR = 'RED'

      ELSE IF (VALUE.EQ.2.) THEN
         CLR = 'GREEN'

      ELSE IF (VALUE.EQ.3.) THEN
         CLR = 'YELLOW'

      ELSE IF (VALUE.EQ.4.) THEN
         CLR = 'BLUE'

      ELSE IF (VALUE.EQ.6.) THEN
         CLR = 'CYAN'

      ELSE IF (VALUE.EQ.5.) THEN
         CLR = 'MAGENTA'

      ELSE IF (VALUE.EQ.7.) THEN
         CLR = 'WHITE'

      ELSE IF (VALUE.EQ.8.) THEN
         CLR = 'GRAY'

      ELSE IF (VALUE.EQ.10.) THEN
         CLR = 'DKGRAY'

      ELSE IF (VALUE.EQ.9.) THEN
         CLR = 'LTGRAY'

      ELSE IF (VALUE.EQ.12.) THEN
         CLR = 'LIME'

      ELSE IF (VALUE.EQ.11.) THEN
         CLR = 'PINK'

      ELSE IF (VALUE.EQ.15.) THEN
         CLR = 'ORANGE'

      ELSE IF (VALUE.EQ.14.) THEN
         CLR = 'VIOLET'

      ELSE IF (VALUE.EQ.13.) THEN
         CLR = 'LTBLUE'

      ELSE
         PLTCNM = .FALSE.
         CLR = 'UNKNOWN'
      END IF

      RETURN

      END
