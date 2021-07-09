C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      INTEGER FUNCTION GETPRC()
C returns the precision for output.
C   GETPRC = 4.  std precision
C   GETPRC = 9.  extended precision
      COMMON /GRPPRC/ IPREC

      GETPRC = IPREC
      RETURN
      END
