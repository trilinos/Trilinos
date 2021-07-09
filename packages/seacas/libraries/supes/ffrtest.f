C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details
      PROGRAM FFRTEST
      PARAMETER (MFIELD=5)
      CHARACTER*16 CVALUE(MFIELD)
      INTEGER KVALUE(MFIELD),IVALUE(MFIELD)
      REAL RVALUE(MFIELD)
      CHARACTER*32 STRING

      CALL GSUPEV(STRING)
      WRITE (*,'(A, A)') 'SUPES VERSION ', STRING
   10 CALL FREFLD( 0,0,'AUTO',MFIELD,IOSTAT,NFIELD,KVALUE,
     *                   CVALUE,IVALUE,RVALUE )
      IF ( IOSTAT .EQ. 0 ) THEN
         PRINT 1000,NFIELD
         PRINT 2000,
     *         (I,KVALUE(I),CVALUE(I),RVALUE(I),IVALUE(I),I=1,MFIELD)
         GO TO 10
      END IF
 1000 FORMAT( ' NFIELD =',I5 /
     *        4X,'I',5X,'KV(I)',15X,'CV(I)',10X,'RV(I)',7X,'IV(I)' )
 2000 FORMAT( I5,I10,'  "',A,'"',1PG15.3,I12 )
      END
