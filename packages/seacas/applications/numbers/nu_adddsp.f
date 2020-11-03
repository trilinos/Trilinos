C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE ADDDSP (COORDS, DSP)
      DIMENSION COORDS (NUMNP,*), DSP(NUMNP,*)
CC
      include 'nu_numg.blk'
      IF (NDIM .EQ. 2) THEN
         DO 10 J=1,NUMNP
            DSP (J, 1) = COORDS (J, 1) + DSP (J, 1)
            DSP (J, 2) = COORDS (J, 2) + DSP (J, 2)
   10    CONTINUE
      ELSE
         DO 20 J=1,NUMNP
            DSP (J, 1) = COORDS (J, 1) + DSP (J, 1)
            DSP (J, 2) = COORDS (J, 2) + DSP (J, 2)
            DSP (J, 3) = COORDS (J, 3) + DSP (J, 3)
   20    CONTINUE
      END IF
CC
      RETURN
      END
