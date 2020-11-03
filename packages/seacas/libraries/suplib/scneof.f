C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SCNEOF
C=======================================================================

C   --*** SCNEOF *** (ETCLIB) Scan input until end of file
C   --   Written by Amy Gilkey - revised 02/23/88
C   --
C   --SCNEOF scans the input file until it reaches the end of file.
C   --Returns without reading if not in batch mode.

      LOGICAL BATCH

      IF (.NOT. BATCH ()) RETURN

      RETURN
      END
