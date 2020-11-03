C    Copyright(C) 1999-2020 National Technology & Engineering Solutions
C    of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C    NTESS, the U.S. Government retains certain rights in this software.
C
C    See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE OPNLOG (LOGU)
C=======================================================================
C   --*** OPNLOG *** (BLOT) Open log file and write header
C   --   Written by Amy Gilkey - revised 12/21/87
C   --
C   --OPNLOG opens the log file and writes the command line as the header
C   --for the log file.
C   --
C   --Parameters:
C   --   NLOG - IN - the log file number
C   --
C   --Common Variables:
C   --   Uses QAINFO of /PROGQA/
C   --   Uses NDBIN, NDBOUT of /DBASE/

      include 'exodusII.inc'
      include 'ag_progqa.blk'
      include 'ag_dbase.blk'

      CHARACTER*256 INLINE
      CHARACTER*256 STR

      NLOG = LOGU
      CALL OPNFIL (NLOG, 'U', 'L', 0, IERR)
      IF (IERR .NE. 0) THEN
         CALL PRTERR ('WARNING', 'Log file cannot be opened')
         NLOG = -1
         GOTO 100
      END IF

      INLINE = '$$$ ' // QAINFO(1)
      L = LENSTR (INLINE) + 1

      CALL EXNAME(NDBIN, STR, LFIL)
      IF (L .LT. LEN (INLINE)) INLINE(L+1:) = STR(:LFIL)
      L = LENSTR (INLINE) + 1

      CALL EXNAME(NDBOUT, STR, LFIL)
      IF (L .LT. LEN (INLINE)) INLINE(L+1:) = STR(:LFIL)
      L = LENSTR (INLINE) + 1

      WRITE (NLOG, '(A)') INLINE(:L-1)

  100 CONTINUE
      LOGU = NLOG
      RETURN
      END
