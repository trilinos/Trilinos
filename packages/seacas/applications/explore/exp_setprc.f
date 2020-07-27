C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      SUBROUTINE SETPRC(I,MSG)
C sets the precision for output.
C   I = 4.  std precision
C   I = 9.  extended precision
C MSG=0 no message
      COMMON /GRPPRC/ IPREC

      IPREC = I
      IF (MSG.NE.0) WRITE(*,10) I
 10   FORMAT('Precision set to ', I2)
      RETURN
      END
