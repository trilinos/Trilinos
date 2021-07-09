C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE OPNLOG (LOGU)
C=======================================================================

C   --*** OPNLOG *** (BLOT) Open log file and write header
C   --   Written by Amy Gilkey - revised 01/11/88
C   --
C   --OPNLOG opens the log file and writes the command line as the header
C   --for the log file.
C   --
C   --Parameters:
C   --   NLOG - IN/OUT - the log file number; returned <= if log file
C   --      cannot be opened
C   --
C   --Common Variables:
C   --   Uses NDB of /DBASE/
C   --   Uses QAINFO of /PROGQA/

      include 'params.blk'
      include 'progqa.blk'
      include 'dbase.blk'
      include 'dbname.blk'

      CHARACTER*2048 INLINE, FILNAM, ERRMSG
      CHARACTER*256 STR
      LOGICAL ISON

      NLOG = LOGU
      filnam = basenam(:lenstr(basenam)) // '.blot.log'

      open (unit=nlog, file=filnam(:lenstr(filnam)), form='formatted',
     *  status='unknown', iostat=ierr)
      IF (IERR .NE. 0) THEN
        ERRMSG = 'Log file "'//FILNAM(:LENSTR(FILNAM))//
     *    '" could not be opened.'
        CALL PRTERR ('CMDERR', ERRMSG(:LENSTR(ERRMSG)))
        GOTO 100
      END IF

      INLINE = '$$$ ' // QAINFO(1)
      L = LENSTR (INLINE) + 1

      IF (L .LT. LEN (INLINE)) INLINE(L+1:) = DBNAME
      L = LENSTR (INLINE) + 1

      CALL GRGPARD ('DEVICE', 1, ISON, STR)
      IF (ISON) THEN
         IF (L .LT. LEN (INLINE)) INLINE(L+1:) = STR
      ELSE
         IF (L .LT. LEN (INLINE)) INLINE(L+1:) = '""'
      END IF
      L = LENSTR (INLINE) + 1

      CALL GRGPARD ('DEVICE', 2, ISON, STR)
      IF (ISON) THEN
         IF (L .LT. LEN (INLINE)) INLINE(L+1:) = STR
      ELSE
         IF (L .LT. LEN (INLINE)) INLINE(L+1:) = '""'
      END IF
      L = LENSTR (INLINE) + 1

      WRITE (NLOG, '(A)', IOSTAT=IERR) INLINE(:L-1)
      if (ierr .ne. 0) then
         CALL PRTERR ('WARNING', 'Log file cannot be written')
         NLOG = -1
         GOTO 100
      end if
  100 CONTINUE
      LOGU = NLOG
      RETURN
      END
