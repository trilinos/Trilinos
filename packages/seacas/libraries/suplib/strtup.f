C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.

C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE STRTUP (QAINFO)
C=======================================================================

C   --*** STRTUP *** (ETCLIB) Startup program
C   --   Written by Amy Gilkey - revised 11/24/87
C   --
C   --STRTUP should be called at the start of any program.  It performs
C   --initialization details common to all programs.  Specifically, it
C   --   Initializes the CPU time
C   --   Sets the program name, time, and date
C   --
C   --Parameters:
C   --   QAINFO - IN/OUT - the current program version information:
C   --      QAINFO(1..3) input, QAINFO(4..6) set
C   --      (1) = program name
C   --      (2) = revision date
C   --      (3) = version as "QA xx.xx" or "X  xx.xx" or "   xx.xx"
C   --      (4) = program name with version appended
C   --      (5) = date of current run
C   --      (6) = time of current run

C   --Routines Called:
C   --   LENSTR - (STRLIB) Find string length

      CHARACTER*(*) QAINFO(6)
      CHARACTER*8   MEMDBG

C   --Initialize the CPU time
      CALL EXCPUS (CPUSEC)

C   --Set the program information array
      QAINFO(4) = QAINFO(1)
      L = MIN (LENSTR(QAINFO(4))+1, 7)
      QAINFO(4)(L:L+1) = QAINFO(3)(4:5)
      IF (QAINFO(4)(L:L) .EQ. ' ') QAINFO(4)(L:L) = '0'
      QAINFO(5) = ' '
      QAINFO(6) = ' '
      CALL EXDATE (QAINFO(5))
      CALL EXTIME (QAINFO(6))

C ... If EXT99 Environment variable set, turn on supes memory debugging
C     The numeric value of the variable is used as the unit to write
C     debug information to.
      CALL EXNAME (-99, MEMDBG, L)
      IF (L .GE. 1) THEN
        READ(MEMDBG(:L), '(I8)', ERR=20) IUNIT
        CALL MDDEBG(IUNIT)
      END IF
 20   CONTINUE

      RETURN
      END
