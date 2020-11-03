C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

      PROGRAM MEMTEST

C     THIS PROGRAM TESTS THE SUPES MEMORY MANAGER.

      PARAMETER (MFIELD=4)
      CHARACTER*32 VERSN
      CHARACTER*8 CV(MFIELD)
      PARAMETER (MAXSIZ=256)
      CHARACTER*1 CBASE(MAXSIZ)
      COMMON /EXTCLB/ CBASE
      DIMENSION KV(MFIELD),IV(MFIELD),RV(MFIELD),BASE(1),IBASE(1)
      EQUIVALENCE (BASE,IBASE)

      CALL GSUPEV(VERSN)
      WRITE (*,'(A, A)') ' SUPES Version ', VERSN

      CALL EXNAME (6, CV(1), LEN)
c      OPEN (6, STATUS='unknown',FILE=CV(1)(1:LEN))
      CALL EXNAME (8, CV(1), LEN)
c      OPEN (8, STATUS='unknown',FILE=CV(1)(1:LEN))
  100 CALL FREFLD( 0,8,'FUNC: ',MFIELD,IOSTAT,N,KV,CV,IV,RV )
      IF ( IOSTAT .NE. 0 ) CV(1) = 'EXIT'
      IF ( CV(1) .EQ. 'MDINIT' ) THEN
         IOFF = IV(2) + 1
         CALL MDINIT (IBASE(IOFF))
      ELSE IF ( CV(1) .EQ. 'MCINIT' ) THEN
         IOFF = IV(2) + 1
         CALL MCINIT (CBASE(IOFF))
      ELSE IF (CV(1) .EQ. 'MDCOMP') THEN
         CALL MDCOMP()
      ELSE IF (CV(1) .EQ. 'MCCOMP') THEN
         CALL MCCOMP()
      ELSE IF (CV(1) .EQ. 'MDDEBG') THEN
         CALL MDDEBG (6)
      ELSE IF (CV(1) .EQ. 'MCDEBG') THEN
         CALL MCDEBG (6)
      ELSE IF (CV(1) .EQ. 'MDDEL') THEN
         CALL MDDEL (CV(2))
      ELSE IF (CV(1) .EQ. 'MCDEL') THEN
         CALL MCDEL (CV(2))
      ELSE IF (CV(1) .EQ. 'MDEFIX') THEN
         CALL MDEFIX (IV(2), IV(3))
         WRITE (6,10000) 'ERROR CODE ', IV(2)
         WRITE (6,10000) 'NEW COUNT  ', IV(3)
         WRITE (8,10000) 'ERROR CODE ', IV(2)
         WRITE (8,10000) 'NEW COUNT  ', IV(3)
      ELSE IF (CV(1) .EQ. 'MCEFIX') THEN
         CALL MCEFIX (IV(2), IV(3))
         WRITE (6,10000) 'ERROR CODE ', IV(2)
         WRITE (6,10000) 'NEW COUNT  ', IV(3)
         WRITE (8,10000) 'ERROR CODE ', IV(2)
         WRITE (8,10000) 'NEW COUNT  ', IV(3)
      ELSE IF (CV(1) .EQ. 'MDEROR') THEN
         CALL MDEROR (6)
         CALL MDEROR (8)
      ELSE IF (CV(1) .EQ. 'MCEROR') THEN
         CALL MCEROR (6)
         CALL MCEROR (8)
      ELSE IF (CV(1) .EQ. 'MDERPT') THEN
         CALL MDERPT (IV(2), IV(3))
         WRITE (6,10000) 'ERROR CODE  ', IV(2)
         WRITE (6,10000) 'ERROR COUNT ', IV(3)
         WRITE (8,10000) 'ERROR CODE  ', IV(2)
         WRITE (8,10000) 'ERROR COUNT ', IV(3)
      ELSE IF (CV(1) .EQ. 'MCERPT') THEN
         CALL MCERPT (IV(2), IV(3))
         WRITE (6,10000) 'ERROR CODE  ', IV(2)
         WRITE (6,10000) 'ERROR COUNT ', IV(3)
         WRITE (8,10000) 'ERROR CODE  ', IV(2)
         WRITE (8,10000) 'ERROR COUNT ', IV(3)
      ELSE IF (CV(1) .EQ. 'MDEXEC') THEN
         WRITE (6,10000) 'POINTER BEFORE ', IP
         WRITE (8,10000) 'POINTER BEFORE ', IP
         CALL MDEXEC ()
         WRITE (6,10000) 'POINTER AFTER ', IP
         WRITE (8,10000) 'POINTER AFTER ', IP
      ELSE IF (CV(1) .EQ. 'MCEXEC') THEN
         WRITE (6,10000) 'POINTER BEFORE ', IP
         WRITE (8,10000) 'POINTER BEFORE ', IP
         CALL MCEXEC ()
         WRITE (6,10000) 'POINTER AFTER ', IP
         WRITE (8,10000) 'POINTER AFTER ', IP
      ELSE IF (CV(1) .EQ. 'MDFIND') THEN
         CALL MDFIND (CV(2), IP, IL)
         WRITE (6, 10000) 'POINTER:',IP,'LENGTH:',IL
         WRITE (8, 10000) 'POINTER:',IP,'LENGTH:',IL
      ELSE IF (CV(1) .EQ. 'MCFIND') THEN
         CALL MCFIND (CV(2), IP, IL)
         WRITE (6, 10000) 'POINTER:',IP,'LENGTH:',IL
         WRITE (8, 10000) 'POINTER:',IP,'LENGTH:',IL
      ELSE IF (CV(1) .EQ. 'MDFILL') THEN
         IF (KV(2) .EQ. 1) THEN
            CALL MDFILL (RV(2))
         ELSE
            TMP = IV(2)
            CALL MDFILL (TMP)
         END IF
      ELSE IF (CV(1) .EQ. 'MCFILL') THEN
         CALL MCFILL (CV(2))
      ELSE IF (CV(1) .EQ. 'MDFOFF') THEN
         CALL MDFOFF ()
      ELSE IF (CV(1) .EQ. 'MCFOFF') THEN
         CALL MCFOFF ()
      ELSE IF (CV(1) .EQ. 'MDGET') THEN
         CALL MDGET (IV(2))
      ELSE IF (CV(1) .EQ. 'MCGET') THEN
         CALL MCGET (IV(2))
      ELSE IF (CV(1) .EQ. 'MDGIVE') THEN
         CALL MDGIVE()
      ELSE IF (CV(1) .EQ. 'MCGIVE') THEN
         CALL MCGIVE()
      ELSE IF (CV(1) .EQ. 'MDLAST') THEN
         CALL MDLAST (IP)
         WRITE (6,10000) 'ERROR NUMBER:',IP
         WRITE (8,10000) 'ERROR NUMBER:',IP
      ELSE IF (CV(1) .EQ. 'MCLAST') THEN
         CALL MCLAST (IP)
         WRITE (6,10000) 'ERROR NUMBER:',IP
         WRITE (8,10000) 'ERROR NUMBER:',IP
      ELSE IF (CV(1) .EQ. 'MDLIST') THEN
         CALL MDLIST (6)
         CALL MDLIST (8)
      ELSE IF (CV(1) .EQ. 'MCLIST') THEN
         CALL MCLIST (6)
         CALL MCLIST (8)
      ELSE IF (CV(1) .EQ. 'MDLONG') THEN
         CALL MDLONG (CV(2), IP, IV(3))
         WRITE (6,10000) 'POINTER:',IP
         WRITE (8,10000) 'POINTER:',IP
      ELSE IF (CV(1) .EQ. 'MCLONG') THEN
         CALL MCLONG (CV(2), IP, IV(3))
         WRITE (6,10000) 'POINTER:',IP
         WRITE (8,10000) 'POINTER:',IP
      ELSE IF (CV(1) .EQ. 'MDMEMS') THEN
         CALL MDMEMS (NSUA, NSUD, NSUV, NSULV)
         WRITE (6,10000) '   ALLOCATED:', NSUA
         WRITE (6,10000) '    DEFERRED:', NSUD
         WRITE (6,10000) '       VOIDS:', NSUV
         WRITE (6,10000) 'LARGEST VOID:', NSULV
         WRITE (8,10000) '   ALLOCATED:', NSUA
         WRITE (8,10000) '    DEFERRED:', NSUD
         WRITE (8,10000) '       VOIDS:', NSUV
         WRITE (8,10000) 'LARGEST VOID:', NSULV
      ELSE IF (CV(1) .EQ. 'MCMEMS') THEN
         CALL MCMEMS (NSUA, NSUD, NSUV, NSULV)
         WRITE (6,10000) '   ALLOCATED:', NSUA
         WRITE (6,10000) '    DEFERRED:', NSUD
         WRITE (6,10000) '       VOIDS:', NSUV
         WRITE (6,10000) 'LARGEST VOID:', NSULV
         WRITE (8,10000) '   ALLOCATED:', NSUA
         WRITE (8,10000) '    DEFERRED:', NSUD
         WRITE (8,10000) '       VOIDS:', NSUV
         WRITE (8,10000) 'LARGEST VOID:', NSULV
      ELSE IF (CV(1) .EQ. 'MDNAME') THEN
         CALL MDNAME (CV(2), CV(3))
      ELSE IF (CV(1) .EQ. 'MCNAME') THEN
         CALL MCNAME (CV(2), CV(3))
      ELSE IF (CV(1) .EQ. 'MDPRNT') THEN
         CALL MDPRNT (CV(2), 6, CV(3))
         CALL MDPRNT (CV(2), 8, CV(3))
      ELSE IF (CV(1) .EQ. 'MCPRNT') THEN
         CALL MCPRNT (CV(2), 6, IV(3))
         CALL MCPRNT (CV(2), 8, IV(3))
      ELSE IF (CV(1) .EQ. 'MDRSRV') THEN
         CALL MDRSRV (CV(2),IP,IV(3))
         WRITE (6,10000) 'POINTER:',IP
         WRITE (8,10000) 'POINTER:',IP
      ELSE IF (CV(1) .EQ. 'MCRSRV') THEN
         CALL MCRSRV (CV(2),IP,IV(3))
         WRITE (6,10000) 'POINTER:',IP
         WRITE (8,10000) 'POINTER:',IP
      ELSE IF (CV(1) .EQ. 'MDSTAT') THEN
         CALL MDSTAT (MNERRS, MNUSED)
         WRITE (6, 10000) 'ERRORS:',MNERRS,'USED:', MNUSED
         WRITE (8, 10000) 'ERRORS:',MNERRS,'USED:', MNUSED
      ELSE IF (CV(1) .EQ. 'MCSTAT') THEN
         CALL MCSTAT (MNERRS, MNUSED)
         WRITE (6, 10000) 'ERRORS:',MNERRS,'USED:', MNUSED
         WRITE (8, 10000) 'ERRORS:',MNERRS,'USED:', MNUSED
      ELSE IF (CV(1) .EQ. 'MDWAIT') THEN
         CALL MDWAIT ()
      ELSE IF (CV(1) .EQ. 'MCWAIT') THEN
         CALL MCWAIT ()
      ELSE IF (CV(1) .EQ. 'MEMY') THEN
         CALL EXMEMY (1, L, M)
      ELSE IF (CV(1) .EQ. 'EXIT') THEN
         CLOSE (6)
         CLOSE (8)
         STOP ' '
      ELSE
         WRITE (6,10000) 'UNKNOWN COMMAND'
         WRITE (8,10000) 'UNKNOWN COMMAND'
      END IF
      GO TO 100
10000 FORMAT ((1X,A,I32))
      END
