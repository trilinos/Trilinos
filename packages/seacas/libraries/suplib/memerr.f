C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MEMERR
C=======================================================================

C   --*** MEMERR *** (ETCLIB) Flag dynamic memory error
C   --   Written by Amy Gilkey - revised 02/23/88
C   --
C   --MEMERR prints an error message for a dynamic memory error.

      CALL MDSTAT (NERR, MEM)
      CALL MDERPT (2, NOVER)
      IF (NERR .LE. NOVER) THEN
         CALL PRTERR ('FATAL', 'Too much dynamic memory requested')
      ELSE
         CALL PRTERR ('PROGRAM', 'Dynamic allocation problem')
         CALL MDEROR (6)
      END IF

      RETURN
      END
