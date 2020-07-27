C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      PROGRAM TSTEXT
      CHARACTER*32 VERSN, LINE
      CHARACTER*8 DATE,TIME,NAME,HARD,SOFT
      character*1 b(5)
      DIMENSION A(5)
      CALL EXCPUS( CPU0 )

      CALL GSUPEV(VERSN)
      WRITE (*,'(A, A)') ' SUPES Version ', VERSN

C ... Commented out since I didn't want to make CMake Test that took input from stdin
C      CALL EXREAD( 'TST: ',LINE,IOSTAT )

      LINE="Gregory Sjaardema"
      CALL EXUPCS( LINE )
      PRINT *,'Input line = ',LINE

      CALL EXDATE( DATE )
      PRINT *,'Date = ',DATE

      CALL EXTIME( TIME )
      PRINT *,'Time = ',TIME

      CALL EXNAME( 1,NAME,LN )
      IF ( LN .EQ. 0 ) THEN
         NAME = ' '
         LN = 1
      END IF
      PRINT *,'Unit 1 name = ',NAME(1:LN)

      CALL EXNAME( 10,NAME,LN )
      IF ( LN .EQ. 0 ) THEN
         NAME = ' '
         LN = 1
      END IF
      PRINT *,'Unit 10 name = ',NAME(1:LN)

      CALL EXNAME( -1,NAME,LN )
      IF ( LN .EQ. 0 ) THEN
         NAME = ' '
         LN = 1
      END IF
      PRINT *,'Symbol 1 = ',NAME(1:LN)

      CALL EXPARM( HARD,SOFT,MODE,KSCU,KNSU,IDAU )
      PRINT *,'Processor = ',HARD,'  System = ',SOFT,'  Mode = ',MODE
      PRINT *,'Character, Numeric, D/A Units: ',KSCU,KNSU,IDAU

      CALL EXMEMY( 10,LOCBLK,MEMRTN )
      PRINT *,'Memory block location and length: ',LOCBLK,MEMRTN

      MEMRTN = -998
      CALL EXMEMY( -10,LOCBLK,MEMRTN )
      PRINT *,'Memory block location and length: ',LOCBLK,MEMRTN

      NN = IXLNUM( A(5) ) - IXLNUM( A )
      PRINT *,'Numeric difference = ',NN

      NN = IXLCHR( B(5) ) - IXLCHR( B )
      PRINT *,'Character difference = ',NN

C ... Burn up some time so we can test the excpus function
      ra = 0.0
      do 100 i=1, 10000000
         ra = ra + sqrt(float(i))
 100  continue
      CALL EXCPUS( CPU1 )
      CPUS = CPU1 - CPU0
      PRINT *,'CPU time = ',CPUS
      END
