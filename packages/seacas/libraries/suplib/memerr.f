C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE MEMERR
C=======================================================================
C$Id: memerr.f,v 1.2 2009/03/25 12:46:02 gdsjaar Exp $
C$Log: memerr.f,v $
CRevision 1.2  2009/03/25 12:46:02  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:15:34  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:15:33  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:37  gdsjaar
c Initial revision
c

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
