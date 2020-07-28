C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SHOWFL (TYPE, NUMESS, IDESS, NEESS, IPEESS, name)
C=======================================================================

      DIMENSION IDESS(*), NEESS(*), IPEESS(*)
      CHARACTER*(*) NAME(*)
      CHARACTER*1 TYPE

      IF (TYPE .EQ. 'S') THEN
      WRITE (*, 20)
      DO 10 I=1, NUMESS
          WRITE (*, 30) I, IDESS(I), NEESS(I), IPEESS(I), NAME(I)
   10 CONTINUE
      ELSE IF (TYPE .EQ. 'N') THEN
      WRITE (*, 40)
      DO 15 I=1, NUMESS
          WRITE (*, 50) I, IDESS(I), NEESS(I), NAME(I)
   15 CONTINUE
      END IF

   20 FORMAT (/' Side Sets:           ID     Elements    ',
     $     '    Nodes')
   40 FORMAT (/' Node Sets:           ID        Nodes')
   30 FORMAT (I11,2X,I11,2X,I11,2X,I11,2X,A)
   50 FORMAT (I11,2X,I11,2X,I11,2X,A)
      RETURN
      END
