C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

      SUBROUTINE GAPOUT (DIRCOS, MASSLV, NUMNIQ, NDIM, IFLGM, IFLGS,
     *   GMTHD)
      include 'nu_io.blk'
      DIMENSION DIRCOS(NDIM + 2,*), MASSLV(2,*)
      CHARACTER*8 GMTHD, COSLAB(3), DISLAB(2)
      LOGICAL ISABRT
      DATA COSLAB /'Cosine X','Cosine Y','Cosine Z'/
      DATA DISLAB /' Normal ','Tangent '/
      LENPAG = 50

      DO 30 IO=IOMIN,IOMAX
         DO 20 IPAG = 1, NUMNIQ, LENPAG
            IF (ISABRT()) RETURN
            WRITE (IO, 40) IFLGM, IFLGS, GMTHD, (COSLAB(I),I=1,NDIM),
     *         DISLAB(1), DISLAB(2)

            DO 10 I=IPAG, MIN(IPAG+LENPAG-1,NUMNIQ)
               WRITE (IO, 50) I, (MASSLV(J,I),J=1,2),
     *            (DIRCOS(J,I),J=1,NDIM+2)
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE

   40 FORMAT ('1',
     *   ' Master Flag = ',I6,', Slave Flag = ',I6,', Method: ',A8,//,
     *   '     #   Master   Slave   ',5(A8,3X))
   50 FORMAT (1X,I5,':',2I8,5(1PE11.3))
      RETURN
      END
