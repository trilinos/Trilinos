C Copyright(C) 1999-2020 National Technology & Engineering Solutions
C of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
C NTESS, the U.S. Government retains certain rights in this software.
C 
C See packages/seacas/LICENSE for details

C=======================================================================
      SUBROUTINE CKEXOD (EXODUS, *)
C=======================================================================
C$Id: ckexod.f,v 1.2 2009/03/25 12:46:01 gdsjaar Exp $
C$Log: ckexod.f,v $
CRevision 1.2  2009/03/25 12:46:01  gdsjaar
CAdd copyright and license notice to all files.
C
CRevision 1.1.1.1  1990/08/14 16:12:06  gdsjaar
CTesting
C
c Revision 1.1  90/08/14  16:12:05  gdsjaar
c Initial revision
c
c Revision 1.1  90/08/09  13:39:04  gdsjaar
c Initial revision
c

C   --*** CKEXOD *** (ETCLIB) Check for EXODUS format
C   --   Written by Amy Gilkey - revised 12/23/87
C   --
C   --CKEXOD prints an error message if the database is not in the EXODUS
C   --database format.
C   --
C   --Parameters:
C   --   EXODUS - IN - true iff EXODUS format
C   --   * - return statement if error

      LOGICAL EXODUS

      IF (.NOT. EXODUS) THEN
         CALL PRTERR ('CMDERR',
     &      'Command not allowed on GENESIS database')
         RETURN 1
      END IF

      RETURN
      END
