C $Id: traniq.f,v 1.1 1991/02/21 15:46:03 gdsjaar Exp $
C $Log: traniq.f,v $
C Revision 1.1  1991/02/21 15:46:03  gdsjaar
C Initial revision
C
      SUBROUTINE TRANIQ (LSTSN, MAP, MASSLV, NSEG, IDIM)
      DIMENSION LSTSN(*), MAP(*), MASSLV(IDIM,*)
C
      DO 10 I=1,NSEG
          MASSLV(1,MAP(I)) = LSTSN(I)
   10 CONTINUE
C
      RETURN
      END
