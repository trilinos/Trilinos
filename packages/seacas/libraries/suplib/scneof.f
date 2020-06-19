C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE SCNEOF
C=======================================================================
C$Id: scneof.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: scneof.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:16:17  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:16:16  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:45  gdsjaar
c Initial revision
c

C   --*** SCNEOF *** (ETCLIB) Scan input until end of file
C   --   Written by Amy Gilkey - revised 02/23/88
C   --
C   --SCNEOF scans the input file until it reaches the end of file.
C   --Returns without reading if not in batch mode.

      LOGICAL BATCH

      IF (.NOT. BATCH ()) RETURN

      RETURN
      END
