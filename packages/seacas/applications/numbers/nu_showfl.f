C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE SHOWFL (TYPE, NUMESS, IDESS, NEESS, IPEESS)
      DIMENSION IDESS(*), NEESS(*), IPEESS(*)
      CHARACTER*1 TYPE

      IF (TYPE .EQ. 'S') THEN
      WRITE (*, 20)
      DO 10 I=1, NUMESS
          WRITE (*, 30) I, IDESS(I), NEESS(I), IPEESS(I)
   10 CONTINUE
      ELSE IF (TYPE .EQ. 'N') THEN
      WRITE (*, 40)
      DO 15 I=1, NUMESS
          WRITE (*, 30) I, IDESS(I), NEESS(I)
   15 CONTINUE
      END IF

   20 FORMAT (/' Side Set Flags:'/
     *    '           ID     Elements    Nodes')
   40 FORMAT (/' Node Set Flags:'/
     *    '           ID     Nodes')
   30 FORMAT (I6,':',3(I6,5X))
      RETURN
      END
