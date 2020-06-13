C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C    
C    See packages/seacas/LICENSE for details

C $Id: showfl.f,v 1.1 1991/02/21 15:45:37 gdsjaar Exp $
C $Log: showfl.f,v $
C Revision 1.1  1991/02/21 15:45:37  gdsjaar
C Initial revision
C
      SUBROUTINE SHOWFL (TYPE, NUMESS, IDESS, NEESS, IPEESS)
      DIMENSION IDESS(*), NEESS(*), IPEESS(*)
      CHARACTER*1 TYPE
C
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
C
   20 FORMAT (/' Side Set Flags:'/
     *    '           ID     Elements    Nodes')
   40 FORMAT (/' Node Set Flags:'/
     *    '           ID     Nodes')
   30 FORMAT (I6,':',3(I6,5X))
      RETURN
      END
